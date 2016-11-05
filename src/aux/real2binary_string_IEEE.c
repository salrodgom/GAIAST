#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

char *fmt_binstr (unsigned long n);
int main (int argc, char** argv) {
    if (argc < 2) {
        fprintf (stderr, "error: insufficient input. Usage: %s float\n", argv[0]);
        return 1;
    }
    unsigned int i = *(unsigned int *)&argv[1];
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
