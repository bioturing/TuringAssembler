#ifndef _KMC_SKIPPING_H_
#define _KMC_SKIPPING_H_

#ifdef __cplusplus
extern "C" {
#endif

int KMC_build_kmer_database(int ksize, const char *working_dir, int n_threads,
					int mmem, int n_files, char **files);

int KMC_arg_kmer_count(int argc, char* argv[]);

#ifdef __cplusplus
}
#endif

#endif
