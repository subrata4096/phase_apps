#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>

long
perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
               int cpu, int group_fd, unsigned long flags)
{
   int ret;

   ret = syscall(__NR_perf_event_open, hw_event, pid, cpu,
                  group_fd, flags);
   return ret;
}

long long perfProfileRead(int* fdptr)
{
   int fd = *fdptr;
   long long count;
   //printf("%d\n", fd);   
   read(fd, &count, sizeof(long long));

   printf("Number of HW instructions Used=%lld\n", count);
   close(fd);

   return count;
}

void perfProfileInit(struct perf_event_attr* peptr, int* fdptr)
{
  int fd = *fdptr;
   //printf("%d\n", fd);   

   peptr->type = PERF_TYPE_HARDWARE;
   peptr->size = sizeof(struct perf_event_attr);
   peptr->config = PERF_COUNT_HW_INSTRUCTIONS;
   peptr->disabled = 1;
   peptr->exclude_kernel = 1;
   peptr->exclude_hv = 1;
   peptr->exclude_idle = 1;

   fd = perf_event_open(peptr, 0, -1, -1, 0);
   if (fd == -1) {
      fprintf(stderr, "Error opening leader %llx\n", peptr->config);
      exit(EXIT_FAILURE);
   }
   *fdptr = fd;
}

void perfProfileStart(int* fdptr)
{
   int fd = *fdptr;
   //printf("%d\n", fd);   
   ioctl(fd, PERF_EVENT_IOC_RESET, 0);
   ioctl(fd, PERF_EVENT_IOC_ENABLE, 0);
}

void perfProfileEnd(int* fdptr)
{
   int fd = *fdptr;
   //printf("%d\n", fd);   
   ioctl(fd, PERF_EVENT_IOC_DISABLE, 0);
}

int
main(int argc, char **argv)
{
   struct perf_event_attr pe;
   int fd;

   memset(&pe, 0, sizeof(struct perf_event_attr));
   perfProfileInit(&pe, &fd);
   perfProfileStart(&fd);

   printf("Measuring instruction count for this printf\n");
   printf("Measuring instruction count for this printf\n");
   printf("Measuring instruction count for this printf\n");

   perfProfileEnd(&fd);
   perfProfileRead(&fd);

   return 0;
}
