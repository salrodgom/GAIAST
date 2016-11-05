#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
extern "C"
{
 char Creal2char_IEEE32(float *);
}
char *fmt_binstr (unsigned long n);
char Creal2char_IEEE32 (float r) {
    unsigned int i = *(unsigned int *)&r;
    printf ("%s\n", fmt_binstr (i));
    return 0;
}
char *fmt_binstr (unsigned long n){
    static char s[128 + 1] = {0};
    char *p = s + 128;
    unsigned char i;
    for (i = 0; i < 32; i++) {
        p--;
        *p = (n >> i & 1) ? '1' : '0';
    }
    return p;
}
