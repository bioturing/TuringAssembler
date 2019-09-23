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
			if (g0.edges[id].seq_len <= STRICT_HEAD_LEN * 2)
				continue;
			char *seq;
			decode_seq(&seq, g0.edges[id].seq, g0.edges[id].seq_len);
			fprintf(fo, ">%d\n%s\n", id, seq);
			free(seq);

			int rc_id = g0.edges[id].rc_id;
			decode_seq(&seq, g0.edges[rc_id].seq, g0.edges[rc_id].seq_len);
			fprintf(fo, ">%d\n%s\n", rc_id, seq);
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
	char *buf[2];
	pthread_mutex_t buf_lock[2];
	pthread_mutex_t file_lock[2];
	int buf_pos[2] = {0};
	FILE *f_reduced[2];
	f_reduced[0] = fopen(reduced_path->R1_path, "wb");
	f_reduced[1] = fopen(reduced_path->R2_path, "wb");
	for (int i = 0; i < 2; ++i){
		pthread_mutex_init(&buf_lock[i], NULL);
		pthread_mutex_init(&file_lock[i], NULL);
		buf[i] = (char *) calloc(BUF_LEN, sizeof(char));
	}
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
		for (int j = 0; j < 2; ++j){
			worker_bundles[i].buf_lock[j] = &buf_lock[j];
			worker_bundles[i].buf[j] = buf[j];
			worker_bundles[i].file_lock[j] = &file_lock[j];
			worker_bundles[i].buf_pos[j] = &buf_pos[j];
			worker_bundles[i].f[j] = f_reduced[j];
		}
	}

	pthread_t *producer_threads = calloc(1, sizeof(pthread_t));
	pthread_t *worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	pthread_create(producer_threads, &attr, fastq_producer, producer_bundles);

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, fastq_reducer_iterator,
				worker_bundles + i);


	pthread_join(producer_threads[0], NULL);
	for (int i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	for (int i = 0; i < 2; ++i){
		if (buf_pos[i] > 0)
			fwrite(buf[i], sizeof(char), buf_pos[i], f_reduced[i]);
		fclose(f_reduced[i]);
	}

	for (int i = 0; i < 2; ++i){
		pthread_mutex_destroy(&buf_lock[i]);
		pthread_mutex_destroy(&file_lock[i]);
		free(buf[i]);
	}
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

	pthread_mutex_t **buf_lock = bundle->buf_lock;
	pthread_mutex_t **file_lock = bundle->file_lock;
	char **buf = bundle->buf;
	int **buf_pos = bundle->buf_pos;
	FILE **f = bundle->f;
	if (is_mapped_head || (!is_mapped_head && !is_mapped_mid)){
		for (int i = 0; i < 2; ++i){
			char own_buf[1024];
			char *name, *info, *qual, *path, *seq;
			if (i == 0){
				name = r1->name;
				info = r1->info;
				qual = r1->qual;
				seq = r1->seq;
			} else {
				name = r2->name;
				info = r2->info;
				qual = r2->qual;
				seq = r2->seq;
			}
			sprintf(own_buf, "@%s\t%s\n%s\n+\n%s\n", name, info, seq,
					qual);
			int len = strlen(own_buf);
			pthread_mutex_lock(buf_lock[i]);
			if (*(buf_pos[i]) + len > BUF_LEN){
				pthread_mutex_lock(file_lock[i]);
				fwrite(buf[i], sizeof(char), *(buf_pos[i]),
						f[i]);
				pthread_mutex_unlock(file_lock[i]);
				*(buf_pos[i]) = 0;
			}
			memcpy(buf[i] + *(buf_pos[i]), own_buf,
					len * sizeof(char));
			*(buf_pos[i]) += len;
			pthread_mutex_unlock(buf_lock[i]);
		}
	}
	free(ar1.a);
	free(ar2.a);
	free(r1_seq);
	free(r2_seq);

}

