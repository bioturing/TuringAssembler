#include "fastq_reducer.h"
#include "barcode_builder.h"
#include "utils.h"
#include "fastq_producer.h"
#include "verbose.h"
#include "assembly_graph.h"
#include "helper.h"

void *fastq_reducer_iterator(void *data)
{
	struct fastq_reducer_bundle_t *bundle = (struct fastq_reducer_bundle_t *)
		data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format, mapper_algo;

	int64_t n_reads;
	int64_t *gcnt_reads;
	uint64_t barcode;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf) {
			break;
		}
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = 0;
		buf1 = ext_buf->R1_buf;
		buf2 = ext_buf->R2_buf;
		input_format = ext_buf->input_format;

		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, buf1, &pos1) :
				get_read_from_fa(&read1, buf1, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, buf2, &pos2) :
				get_read_from_fa(&read2, buf2, &pos2);


			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				__ERROR("\nWrong format file\n");

			reduce_read(&read1, &read2, bundle);

			if (rc1 == READ_END)
				break;
		}
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

void fastq_reducer(struct opt_proc_t *opt, struct read_path_t *org_rpath,
		struct read_path_t *reduced_path)
{
	char fasta_path[1024];
	sprintf(fasta_path, "%s/all_scaffold_contigs.fasta", opt->out_dir);
	struct asm_graph_t g0;
	load_asm_graph(&g0, opt->in_file);
	FILE *f = fopen(opt->in_fasta, "r");
	FILE *fo = fopen(fasta_path, "w");
	int n_paths;
	fscanf(f, "%d\n", &n_paths);
	for (int i = 0; i < n_paths; ++i){
		int n_contigs;
		fscanf(f, "%d\n", &n_contigs);
		for (int j = 0; j < n_contigs; ++j){
			int id;
			fscanf(f, "%d ", &id);
			char *seq;
			decode_seq(&seq, g0.edges[id].seq, g0.edges[id].seq_len);
			fprintf(fo, ">%d\n%s\n", id, seq);
			free(seq);
		}
	}
	fclose(f);
	fclose(fo);

	bwa_idx_build(fasta_path, fasta_path, BWTALGO_AUTO, 500000000);
	bwaidx_t *bwa_idx = bwa_idx_load(fasta_path, BWA_IDX_ALL);
	mem_opt_t *bwa_opt = asm_memopt_init();
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	struct producer_bundle_t *producer_bundles = init_fastq_pair(opt->n_threads,
			1, &(org_rpath->R1_path), &(org_rpath->R2_path));

	struct fastq_reducer_bundle_t *worker_bundles = calloc(opt->n_threads,
			sizeof(struct fastq_reducer_bundle_t));
	for (int i = 0; i < opt->n_threads; ++i){
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].bwa_idx = bwa_idx;
		worker_bundles[i].bwa_opt = bwa_opt;
		worker_bundles[i].out_rpath = reduced_path;
	}

	pthread_t *producer_threads = calloc(1, sizeof(pthread_t));
	pthread_t *worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	pthread_create(producer_threads, &attr, fastq_producer, producer_bundles);

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, fastq_reducer_iterator,
				worker_bundles + i);


	pthread_join(producer_threads[0], NULL);

	// TODO 
	free_fastq_pair(producer_bundles, 1);
	free(worker_bundles);
	free(producer_threads);
	free(worker_threads);
	bwa_idx_destroy(bwa_idx);
	free(bwa_opt);
	asm_graph_destroy(&g0);
}

void reduce_read(struct read_t *r1, struct read_t *r2,
		struct fastq_reducer_bundle_t *bundle)
{
	bwaidx_t *idx = bundle->bwa_idx;
	mem_opt_t *opt = bundle->bwa_opt;
	mem_alnreg_v ar1, ar2;
	uint8_t *r1_seq, *r2_seq;
	int n1, n2, count, best_score_1, best_score_2;
	ar1 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r1->len, r1->seq);
	ar2 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r2->len, r2->seq);
	r1_seq = malloc(r1->len);
	r2_seq = malloc(r2->len);
	for (int i = 0; i < r1->len; ++i)
		r1_seq[i] = nst_nt4_table[(int)r1->seq[i]];
	for (int i = 0; i < r2->len; ++i)
		r2_seq[i] = nst_nt4_table[(int)r2->seq[i]];
	struct asm_align_t *p1, *p2;
	p1 = alloca(ar1.n * sizeof(struct asm_align_t));
	p2 = alloca(ar2.n * sizeof(struct asm_align_t));
	n1 = n2 = 0;
	best_score_1 = best_score_2 = - (1 << 30);
	for (int i = 0; i < (int)ar1.n; ++i) {
		struct asm_align_t a;
		a = asm_reg2aln(opt, idx->bns, idx->pac, r1->len, r1_seq, ar1.a + i);
		if (a.rid == -1)
			continue;
		if (a.score > best_score_1) {
			best_score_1 = a.score;
			p1[0] = a;
			n1 = 1;
		} else if (a.score == best_score_1) {
			p1[n1++] = a;
		}
	}
	for (int i = 0; i < (int)ar2.n; ++i) {
		struct asm_align_t a;
		a = asm_reg2aln(opt, idx->bns, idx->pac, r2->len, r2_seq, ar2.a + i);
		if (a.rid == -1)
			continue;
		if (a.score > best_score_2) {
			best_score_2 = a.score;
			p2[0] = a;
			n2 = 1;
		} else if (a.score == best_score_2) {
			p2[n2++] = a;
		}
	}
	int is_mapped_mid = 0;
	int is_mapped_head = 0;
	for (int i = 0; i < n1; ++i) {
		if (p1[i].aligned < r1->len)
			continue;
		int pos = p1[i].pos;
		int ref_len = idx->bns->anns[p1[i].rid].len;
		if (pos <= STRICT_HEAD_LEN || pos >= ref_len - STRICT_HEAD_LEN){
			is_mapped_head = 1;
			break;
		} else {
			is_mapped_mid = 1;
		}
	}
	for (int i = 0; i < n2; ++i){
		if (p2[i].aligned < r2->len)
			continue;
		int pos = p2[i].pos;
		int ref_len = idx->bns->anns[p2[i].rid].len;
		if (pos <= STRICT_HEAD_LEN || pos >= ref_len - STRICT_HEAD_LEN){
			is_mapped_head = 1;
			break;
		} else {
			is_mapped_mid = 1;
		}
	}

	if (is_mapped_head || (!is_mapped_head && !is_mapped_mid)){
		__VERBOSE("%s\n%s\n", r1_seq, r2_seq);
	}
	free(ar1.a);
	free(ar2.a);
	free(r1_seq);
	free(r2_seq);

}

