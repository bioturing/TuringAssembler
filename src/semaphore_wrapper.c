#include <stdint.h>
#include "semaphore_wrapper.h"

inline void sem_wrap_init(struct sem_wrap_t *sem, uint32_t value)
{
#if defined (__APPLE__) && defined (__MACH__)
	sem->sem = dispatch_semaphore_create(value);
#else
	sem_init(&sem->sem, 0, value);
#endif
}

inline void sem_wrap_wait(struct sem_wrap_t *sem)
{
#if defined (__APPLE__) && defined (__MACH__)
	dispatch_semaphore_wait(sem->sem, DISPATCH_TIME_FOREVER);
#else
	sem_wait(&sem->sem);
#endif
}

inline void sem_wrap_post(struct sem_wrap_t *sem)
{
#if defined (__APPLE__) && defined (__MACH__)
	dispatch_semaphore_signal(sem->sem);
#else
	sem_post(&sem->sem);
#endif
}

inline void sem_wrap_destroy(struct sem_wrap_t *sem)
{
#if defined (__APPLE__) && defined (__MACH__)
	dispatch_release(sem->sem);
#else
	sem_destroy(&sem->sem);
#endif
}

