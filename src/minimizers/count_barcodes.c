/*
 * =====================================================================================
 *
 *       Filename:  count_barcodes.c
 *
 *    Description:  Count the frequency of each barcode
 *
 *        Version:  1.0
 *        Created:  09/22/2019 10:55:15 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Phan Thai Duy Tan (ptdtan), tan@bioturing.com
 *   Organization:  BioTuring Dev Team
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fastq_producer.h"

KHASH_MAP_INIT_INT(khBx, int);
static khash_t(khBx) bx_freq;

void count_bx_freq(struct opt_proc_t *opt, struct read_path_t *r_path)
{
			
}

static void *biot_buffer_iterator_simple(void *data)
{
	struct readsort_bundle_t *bundle = (struct readsort_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2, readbc;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_trip_buffer();
	char *prefix = bundle->prefix;
	int64_t sm = bundle->sm;

	char *buf, *fbuf;
	int64_t buf_len, buf_size;
	buf = malloc(sm);
	buf_len = 0;
	buf_size = sm;

	char *R1_buf, *R2_buf;
	int pos1, pos2, rc1, rc2, input_format;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf)
			break;
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = 0;
		R1_buf = ext_buf->R1_buf;
		R2_buf = ext_buf->R2_buf;
		input_format = ext_buf->input_format;

		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, R1_buf, &pos1) :
				get_read_from_fa(&read1, R1_buf, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, R2_buf, &pos2) :
				get_read_from_fa(&read2, R2_buf, &pos2);

			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				__ERROR("\nWrong format file\n");

			/* read_name + \t + BX:Z: + barcode + \t + QB:Z: + barcode_quality + \n */
			uint64_t barcode = get_barcode_biot(read1.info, &readbc);
			int record_len, len1, len2;
			if (barcode != (uint64_t)-1) {
			} else {
			}
			if (rc1 == READ_END)
				break;
		}
	}
	return NULL;
}
