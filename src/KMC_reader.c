#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "attribute.h"
#include "assembly_graph.h"
#include "io_utils.h"
#include "KMC_reader.h"
#include "utils.h"
#include "verbose.h"
#include "log.h"

void destroy_kmc_info(struct kmc_info_t *inf)
{
	free(inf->prefix_file_buf);
	free(inf->signature_map);
	inf->prefix_file_buf = NULL;
	inf->signature_map = NULL;
}

typedef unsigned long long uint64;
void KMC_read_prefix(const char *path, struct kmc_info_t *data)
{
	FILE *fp = xfopen(path, "rb");
	char sig[4];
	/* read signature to check whether not truncated file */
	fseek(fp, 0, SEEK_SET);
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, KMCP, 4))
		__ERROR("Wrong KMC prefix file");
	fseek(fp, -4, SEEK_END);
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, KMCP, 4))
		__ERROR("KMC prefix file is truncated");
	fseek(fp, 0, SEEK_END);
	uint64_t size = ftell(fp)-4-4; // size of file without 2 marker and header_offset

	/* read header position */
	int header_offset;
	fseek(fp, -8, SEEK_END);
	xfread(&header_offset, 4, 1, fp);
	log_debug("header position %d", header_offset);

	int kmer_type;
	fseek(fp, -12, SEEK_END);
	xfread(&kmer_type, 4, 1, fp);
	fseek(fp, ftell(fp), SEEK_SET);
	fseek(fp, 4, SEEK_SET);

	if (kmer_type == 0x200) {
	    size -= 4;
        /* read header */
        fseek(fp, -header_offset - 8, SEEK_END);
        xfread(&(data->header), sizeof(struct kmc_header_t), 1, fp);

        /* calculate some auxiliary information */
        data->signature_map_size = ((uint64_t) 1 << (2 * data->header.signature_length)) + 1;
        uint64_t sig_map_size_byte = data->signature_map_size * sizeof(uint32_t);
        uint64_t lut_area_size_in_bytes = size - sig_map_size_byte - header_offset - 8;
        data->prefix_file_buf_size = (lut_area_size_in_bytes + 8) / sizeof(uint64_t);

        /* read prefix offset */
        data->prefix_file_buf= malloc(data->prefix_file_buf_size*8);
        fseek(fp, 4, SEEK_SET);
        size_t result = xfread(data->prefix_file_buf, 1, lut_area_size_in_bytes+8, fp);
        if (result == 0)
            log_error("wtf am i doing here");

        data->prefix_file_buf[data->prefix_file_buf_size] = data->header.total_kmers + 1;

        /* read signature map */
        data->signature_map = malloc(sig_map_size_byte);
        xfread(data->signature_map, data->signature_map_size, sizeof(uint32_t), fp);
        fclose(fp);
    } else if (kmer_type == 0) {
	    log_debug("SIZE %d", size);
        data->prefix_file_buf_size = (size - 4) / 8;		//reads without 4 bytes of a header_offset (and without markers)
        data->prefix_file_buf = malloc(8*data->prefix_file_buf_size);
        size_t result = xfread(data->prefix_file_buf, 8, data->prefix_file_buf_size, fp);
        log_debug("I can read", size);
        if (result == 0)
            log_error("wtf i am doing here, too");

        fseek(fp, -8, SEEK_END);

        size = size - 4;

        unsigned long long header_index = (size - header_offset) / sizeof(uint64);
        uint64 last_data_index = header_index;

        uint64 d = data->prefix_file_buf[header_index];

        data->header.kmer_length = (uint32_t)d;			//- kmer's length
        log_debug("kmer length load: %d", data->header.kmer_length);
        data->header.mode = d >> 32;				//- mode: 0 or 1
        log_debug("mode: %d", data->header.mode);

        header_index++;
        data->header.counter_size = (uint32_t)data->prefix_file_buf[header_index];	//- the size of a counter in bytes;
        log_debug("counter size: %d", data->header.counter_size);
        //- for mode 0 counter_size is 1, 2, 3, or 4 (or 5, 6, 7, 8 for small k values)
        //- for mode = 1 counter_size is 4;
        data->header.lut_prefix_length = data->prefix_file_buf[header_index] >> 32;		//- the number of prefix's symbols cut frm kmers;
        log_debug("lut prefix len: %d", data->header.lut_prefix_length);
        //- (kmer_length - lut_prefix_length) is divisible by 4

        header_index++;
        data->header.min_count = (uint32_t)data->prefix_file_buf[header_index];    //- the minimal number of kmer's appearances
        log_debug("min count : %d", data->header.min_count);
        data->header.max_count = data->prefix_file_buf[header_index] >> 32;      //- the maximal number of kmer's appearances

        header_index++;
        data->header.total_kmers = data->prefix_file_buf[header_index];					//- the total number of kmers

        header_index++;
        data->header.both_strands = (data->prefix_file_buf[header_index] & 0x000000000000000F) == 1;
        data->header.both_strands = !data->header.both_strands;

        data->header.max_count += data->prefix_file_buf[header_index] & 0xFFFFFFFF00000000;
        log_debug("max count : %d", data->header.max_count);

        data->prefix_file_buf[last_data_index] = data->header.total_kmers + 1;
        log_debug("prove total kmer load: %lld", data->header.total_kmers);

//        data->header.sufix_size = (data->header.kmer_length - data->header.lut_prefix_length) / 4;
//
//        sufix_rec_size = sufix_size + counter_size;
//
//        return true;
//
//        fseek(fp, -header_offset - 8, SEEK_END);
//        xfread(&(data->header.kmer_length), 4, 1, fp);
//        log_debug("kmer length load: %d\n", data->header.kmer_length);
//        xfread(&(data->header.mode), 4, 1, fp);
//        log_debug("mode: %d\n", data->header.mode);
//        xfread(&(data->header.counter_size), 4, 1, fp);
//        log_debug("counter size: %d\n", data->header.counter_size);
//        xfread(&(data->header.lut_prefix_length), 4, 1, fp);
//        log_debug("lut prefix len: %d\n", data->header.lut_prefix_length);
//        xfread(&(data->header.min_count), 4, 1, fp);
//        log_debug("min count : %d\n", data->header.min_count);
//        xfread(&(data->header.max_count), 4, 1, fp);
//        log_debug("max count : %d\n", data->header.max_count);
//        xfread(&(data->header.total_kmers), 8, 1, fp);
//        log_debug("total kmer load: %lld\n", data->header.total_kmers);
//        fclose(fp);
	}else {
	    log_error("wrong format file kmc read prefix");
	}
}

struct reader_buffer_t {
	int buf_size;
	int cur_pos;
	int rem_size;
	uint8_t *buf;
};

void init_reader_buffer(struct reader_buffer_t *b, int size)
{
	b->buf_size = size;
	b->cur_pos = 0;
	b->rem_size = 0;
	b->buf = malloc(size);
}

void destroy_reader_buffer(struct reader_buffer_t *b)
{
	free(b->buf);
	b->buf = NULL;
	b->buf_size = 0;
}

static inline void KMC_add_prefix_kmer(uint8_t *kmer, int k, int len,
								uint64_t prefix)
{
	int l = (len + 3) / 4, i;
	for (i = 0; i < l; ++i) {
		kmer[k++] = prefix & 0xff;
		prefix >>= 8;
	}
}

static inline void fill_kmer_char(uint8_t k, char *b, int len)
{
	while (len) {
		b[--len] = nt4_char[(uint32_t)(k & 0x3)];
		k >>= 2;
	}
}

struct KMC_bundle_t {
	pthread_mutex_t *lock;
	const char *path;
	struct kmc_info_t *data;
	void *bundle;
	void (*process)(int, uint8_t *, uint32_t, void*);
	int thread_no;
	uint64_t lo_prefix;
	uint64_t hi_prefix;
	uint64_t offset_kmer;
};

void *KMC_worker_multi(void *raw_data)
{
	struct KMC_bundle_t *bundle = (struct KMC_bundle_t *)raw_data;
	pthread_mutex_t *lock = bundle->lock;
	struct kmc_info_t *data = bundle->data;
	void (*process)(int, uint8_t *, uint32_t, void *) = bundle->process;
	void *child_bundle = bundle->bundle;
	int thread_no = bundle->thread_no;
	uint64_t lo_prefix = bundle->lo_prefix;
	uint64_t hi_prefix = bundle->hi_prefix;
	int suffix_size = (data->header.kmer_length - data->header.lut_prefix_length) / 4;
	int record_size = suffix_size + data->header.counter_size;

	FILE *fp = xfopen(bundle->path, "rb");
	fseek(fp, 4 + bundle->offset_kmer * record_size, SEEK_SET);
	uint64_t prefix_mask, i, prefix, n_kmers, k;
	prefix_mask = ((uint64_t) 1 << (2 * data->header.lut_prefix_length)) - 1;
	struct reader_buffer_t b;
	init_reader_buffer(&b, 1 << 24);
	b.rem_size = fread(b.buf, 1, b.buf_size, fp);
	b.cur_pos = 0;
	uint8_t *kmer = alloca((data->header.kmer_length + 3) / 4);
	// __VERBOSE("lo_prefix = %lu; hi_prefix = %lu\n", lo_prefix, hi_prefix);
	for (i = lo_prefix; i < hi_prefix; ++i) {
		prefix = i & prefix_mask;
		n_kmers = data->prefix_file_buf[i + 1] - data->prefix_file_buf[i];
		for (k = 0; k < n_kmers; ++k) {
			/* fill buffer */
			if (b.rem_size < record_size) {
				memmove(b.buf, b.buf + b.cur_pos, b.rem_size);
				b.rem_size += fread(b.buf + b.rem_size, 1,
						b.buf_size - b.rem_size, fp);
				b.cur_pos = 0;
			}
			if (b.rem_size < record_size)
				__ERROR("Truncated KMC suffix file");
			int j;
			for (j = 0; j < suffix_size; ++j)
				kmer[suffix_size - j - 1] = b.buf[b.cur_pos + j];
			KMC_add_prefix_kmer(kmer, suffix_size,
					data->header.lut_prefix_length, prefix);
			uint32_t counter = *((uint32_t *)(b.buf + b.cur_pos + suffix_size));
			// pthread_mutex_lock(lock);
			process(thread_no, kmer, counter, child_bundle);
			// pthread_mutex_unlock(lock);
			b.cur_pos += record_size;
			b.rem_size -= record_size;
		}
	}
	destroy_reader_buffer(&b);
	fclose(fp);
	return NULL;
}

void KMC_retrieve_kmer_multi(const char *path, int n_threads,
			struct kmc_info_t *data, void *bundle,
			void (*process)(int, uint8_t *, uint32_t, void*))
{
	FILE *fp = fopen(path, "rb");
	char sig[4];
	/* read signature to check whether truncated file */
	fseek(fp, 0, SEEK_SET);
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, KMCS, 4))
		__ERROR("Wrong KMC suffix format file");
	fseek(fp, -4, SEEK_END);
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, KMCS, 4))
		__ERROR("[Open] Truncated KMC suffix file");

	/* Calculate record size */
	fseek(fp, 0, SEEK_END);
	uint64_t file_size = ftell(fp) - 8;
	if (file_size % data->header.total_kmers != 0)
		__ERROR("Total kmer does not divide kmer record field total size");
	int record_size = file_size / data->header.total_kmers;
	int suffix_size = (data->header.kmer_length - data->header.lut_prefix_length) / 4;
	if (record_size != suffix_size + (int)data->header.counter_size)
		__ERROR("KMC record size not consistent (%d != %d + %u)",
			record_size, suffix_size, data->header.counter_size);
	fclose(fp);

	/* calculate break of each thread */
	uint64_t i, cnt_kmers, cap, cur_prefix, cur_cnt;
	struct KMC_bundle_t *worker_bundle = calloc(n_threads, sizeof(struct KMC_bundle_t));
	cap = data->header.total_kmers / n_threads + 1;
	cur_prefix = 0;
	cur_cnt = 0;
	cnt_kmers = 0;
	int k = 0;
	i = 0;
	worker_bundle[0].lo_prefix = 0;
	worker_bundle[0].offset_kmer = 0;
	while (data->prefix_file_buf[i + 1] != data->header.total_kmers + 1) {
		cnt_kmers += data->prefix_file_buf[i + 1] - data->prefix_file_buf[i];
		if (cnt_kmers >= cap * (k + 1)) {
			// worker_bundle[k].lo_prefix = cur_prefix;
			worker_bundle[k].hi_prefix = i + 1;
			worker_bundle[k + 1].lo_prefix = i + 1;
			worker_bundle[k + 1].offset_kmer = cnt_kmers;
			// cur_prefix = i + 1;
			// cur_cnt = cnt_kmers;
			++k;
		}
		++i;
	}
	worker_bundle[n_threads - 1].hi_prefix = i;
	// pthread_mutex_t lock;
	// pthread_mutex_init(&lock, NULL);
	for (k = 0; k < n_threads; ++k) {
		// worker_bundle[k].lock = &lock;
		worker_bundle[k].thread_no = k;
		worker_bundle[k].path = path;
		worker_bundle[k].data = data;
		worker_bundle[k].bundle = bundle;
		worker_bundle[k].process = process;
	}
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_t *threads = calloc(n_threads, sizeof(pthread_t));
	for (k = 0; k < n_threads; ++k)
		pthread_create(threads + k, &attr, KMC_worker_multi, worker_bundle + k);
	for (k = 0; k < n_threads; ++k)
		pthread_join(threads[k], NULL);
	// pthread_mutex_destroy(&lock);
	free(threads);
	free(worker_bundle);
}

void KMC_retrive_kmer(const char *path, struct kmc_info_t *data, void *bundle,
				void (*process)(uint8_t *, uint32_t, void *))
{
	FILE *fp = fopen(path, "rb");
	char sig[4];
	/* read signature to check whether truncated file */
	fseek(fp, 0, SEEK_SET);
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, KMCS, 4))
		__ERROR("Wrong KMC suffix format file");
	fseek(fp, -4, SEEK_END);
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, KMCS, 4))
		__ERROR("[Open] Truncated KMC suffix file");

	/* Calculate record size */
	fseek(fp, 0, SEEK_END);
	uint64_t file_size = ftell(fp) - 8;
	if (file_size % data->header.total_kmers != 0)
		__ERROR("Total kmer does not divide kmer record field total size");
	int record_size = file_size / data->header.total_kmers;
	int suffix_size = (data->header.kmer_length - data->header.lut_prefix_length) / 4;
	if (record_size != suffix_size + (int)data->header.counter_size)
		__ERROR("KMC record size not consistent (%d != %d + %u)",
			record_size, suffix_size, data->header.counter_size);

	/* iterate through buffer */
	uint64_t prefix, prefix_mask, n_kmers, k, i;
	prefix_mask = ((uint64_t)1 << (2 * data->header.lut_prefix_length)) - 1;
	i = 0;
	struct reader_buffer_t b;
	init_reader_buffer(&b, 1 << 24);
	fseek(fp, 4, SEEK_SET);
	b.rem_size = fread(b.buf, 1, b.buf_size, fp);
	b.cur_pos = 0;
	uint8_t *kmer = alloca((data->header.kmer_length + 3) / 4);
	while (data->prefix_file_buf[i + 1] != data->header.total_kmers + 1) {
		prefix = i & prefix_mask;
		n_kmers = data->prefix_file_buf[i + 1] - data->prefix_file_buf[i];
		for (k = 0; k < n_kmers; ++k) {
			/* fill buffer */
			if (b.rem_size < record_size) {
				memmove(b.buf, b.buf + b.cur_pos, b.rem_size);
				b.rem_size += fread(b.buf + b.rem_size, 1,
						b.buf_size - b.rem_size, fp);
				b.cur_pos = 0;
			}
			if (b.rem_size < record_size)
				__ERROR("Truncated KMC suffix file");
			int j;
			for (j = 0; j < suffix_size; ++j)
				kmer[suffix_size - j - 1] = b.buf[b.cur_pos + j];
			KMC_add_prefix_kmer(kmer, suffix_size,
					data->header.lut_prefix_length, prefix);
			uint32_t counter = *((uint32_t *)(b.buf + b.cur_pos + suffix_size));
			process(kmer, counter, bundle);
			b.cur_pos += record_size;
			b.rem_size -= record_size;
		}
		++i;
	}
	destroy_reader_buffer(&b);
	fclose(fp);
}

