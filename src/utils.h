/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#ifndef UTILS_H
#define UTILS_H

#include <stdint.h>
#include <stdio.h>
#include <zlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>

#ifndef LC_kroundup32
#define LC_kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#ifndef LC_kroundup64
#define LC_kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#endif
 
#ifdef __GNUC__
// Tell GCC to validate printf format string and args
#define ATTRIBUTE(list) __attribute__ (list)
#else
#define ATTRIBUTE(list)
#endif

#define err_fatal_simple(msg) _err_fatal_simple(__func__, msg)
#define err_fatal_simple_core(msg) _err_fatal_simple_core(__func__, msg)

#define _err_func_printf(fmt, ...) err_func_format_printf(__func__, fmt, ##__VA_ARGS__)
#define _err_color_printf(type, fmt, ...) err_color_format_printf(type, fmt, ##__VA_ARGS__)
#define _err_warning(fmt, ...) _err_color_printf('W', fmt, ##__VA_ARGS__)
#define _err_error(fmt, ...) _err_color_printf('E', fmt, ##__VA_ARGS__)
#define _err_error_exit(fmt, ...) { _err_color_printf('E', fmt, ##__VA_ARGS__) ; exit(EXIT_FAILURE); }
#define _err_success(fmt, ...) _err_color_printf('S', fmt, ##__VA_ARGS__)
#define _err_info(fmt, ...) _err_color_printf('I', fmt, ##__VA_ARGS__)
#define _err_cmd(fmt, ...) _err_color_printf('C', fmt, ##__VA_ARGS__)

#define xopen(fn, mode) err_xopen_core(__func__, fn, mode)
#define xreopen(fn, mode, fp) err_xreopen_core(__func__, fn, mode, fp)
#define xzopen(fn, mode) err_xzopen_core(__func__, fn, mode)

#define xassert(cond, msg) if ((cond) == 0) _err_fatal_simple_core(__func__, msg)

#define _err_simple_func_printf(msg) err_func_printf(__func__, msg)

typedef struct {
	uint64_t x, y;
} pair64_t;

typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; pair64_t *a; } pair64_v;

#ifdef __cplusplus
extern "C" {
#endif

	void err_fatal(const char *header, const char *fmt, ...) ATTRIBUTE((noreturn));
	void err_fatal_core(const char *header, const char *fmt, ...) ATTRIBUTE((noreturn));
	void _err_fatal_simple(const char *func, const char *msg) ATTRIBUTE((noreturn));
	void _err_fatal_simple_core(const char *func, const char *msg) ATTRIBUTE((noreturn));
	FILE *err_xopen_core(const char *func, const char *fn, const char *mode);
	FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp);
	gzFile err_xzopen_core(const char *func, const char *fn, const char *mode);
    size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);
	size_t err_fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream);

	int err_gzread(gzFile file, void *ptr, unsigned int len);
	int err_fseek(FILE *stream, long offset, int whence);
#define err_rewind(FP) err_fseek((FP), 0, SEEK_SET)
	long err_ftell(FILE *stream);
	int err_fprintf(FILE *stream, const char *format, ...)
        ATTRIBUTE((format(printf, 2, 3)));
	int err_printf(const char *format, ...)
        ATTRIBUTE((format(printf, 1, 2)));
	int err_func_printf(const char *func, const char *format, ...)
        ATTRIBUTE((format(printf, 2, 3)));
	int stdout_printf(const char *format, ...)
        ATTRIBUTE((format(printf, 1, 2)));
	int err_fputc(int c, FILE *stream);
#define err_putchar(C) err_fputc((C), stdout)
	int err_fputs(const char *s, FILE *stream);
	int err_puts(const char *s);
    void err_fgets(char *buff, size_t s, FILE *fp);
	int err_fflush(FILE *stream);
	int err_fclose(FILE *stream);
	int err_gzclose(gzFile file);

#define _err_malloc(s) err_malloc(__func__, s)
#define _err_calloc(n, s) err_calloc(__func__, n, s)
#define _err_realloc(p, s) err_realloc(__func__, p, s)
    void *err_malloc(const char* func, size_t s);
    void *err_calloc(const char* func, size_t n, size_t s);
    void *err_realloc(const char* func, void *p, size_t s);

    void usr_sys_cputime(double *usr_t, double *sys_t);
	double cputime();
	double realtime();
    long peakrss(void);
    void print_format_time(FILE *out);
    int err_func_format_printf(const char *func, const char *format, ...);
    int err_color_format_printf(const char type, const char *format, ...);

	void ks_introsort_64 (size_t n, uint64_t *a);
	void ks_introsort_128(size_t n, pair64_t *a);

    char *retrieve_full_url(const char *short_url);


#ifdef __cplusplus
}
#endif

#define _uni_realloc(p, n, m, type) {                   \
    if (m <= 0) {                                       \
        m = 1;                                          \
        m = MAX_OF_TWO(n, m);                           \
        p = (type*)_err_malloc((m) * sizeof(type));     \
    } else if (n >= m) {                                \
        m = n + 1; kroundup32(m);                       \
        p = (type*)_err_realloc(p, (m) * sizeof(type)); \
    }                                                   \
}

#define _realloc(p, m, type) {(m) <<= 1; p = (type*)_err_realloc(p, (m) * sizeof(type));}

#define _sim_insert_abpoa_utils(v, p, n, m, type) { \
    if (n == m) {               \
        _realloc(p, m, type)    \
    }                           \
    p[n++] = v;                 \
}

#define _insert_abpoa_utils(v, p, n, m, type) { \
    int _i, _flag=0;                  \
    for (_i = 0; _i < n; ++_i) {       \
        if (p[_i] == v) {            \
            _flag = 1;               \
            break;                  \
        }                           \
    }                               \
    if (_flag == 0) {                \
        if (n == m) {               \
            _realloc(p, m, type)    \
        }                           \
        p[n++] = v;                 \
    }                               \
}

#define _bin_insert_abpoa_utils_idx(v, p, n, m, type, flag, k_i) { \
    flag=0, k_i=-1;   \
    int _left=0,_right=n-1,_mid;    \
    type _mid_v, _tmp_v;                 \
    if (_right == -1) k_i = 0;   \
    else {                      \
        while (_left <= _right) { \
            _mid = (_left+_right) >> 1;    \
            _mid_v = p[_mid];             \
            if (_mid_v == v) {           \
                k_i = _mid; \
                flag = 1; break;        \
            } else if (_mid_v > v) {     \
                if (_mid != 0) {         \
                    _tmp_v = p[_mid-1];   \
                }                       \
                if (_mid == 0 || v > _tmp_v) { \
                    k_i = _mid;          \
                    break;              \
                }                       \
                else _right = _mid-1;     \
            } else _left = _mid+1;        \
        }                               \
    }                                   \
    if (k_i == -1) k_i = n;         \
}
     
#define _bin_insert_abpoa_utils(v, p, n, m, type) { \
    int _k_i, _flag;    \
    _bin_insert_abpoa_utils_idx(v, p, n, m, type, _flag, _k_i)   \
    if (_flag == 0) {                \
        if (n == m) {               \
            _realloc(p, m, type)    \
        }                           \
        if (_k_i <= n-1)             \
            memmove(p+_k_i+1, p+_k_i, (n-_k_i)*sizeof(type));  \
        (p)[_k_i] = v;               \
        (n)++;                      \
    }                               \
}

#define _bin_search(v, p, n, type, hit, i) { \
    int _left =0,_right=n-1,_mid; \
    type _mid_v;    \
    hit = 0;               \
    if (_right == -1) hit=0;   \
    else {  \
        while (_left <= _right) {   \
            _mid = (_left+_right) >> 1; \
            _mid_v = p[_mid];       \
            if (_mid_v == v) {  \
                i = _mid;   \
                hit = 1;   \
                break;      \
            } else if (_mid_v > v) {   \
                _right = _mid-1;    \
            } else {    \
                _left = _mid+1; \
            }   \
        }   \
    }   \
} 

#define MIN_OF_TWO(a, b) ((a) < (b) ? (a) : (b))
#define MAX_OF_TWO(a, b) ((a) > (b) ? (a) : (b))
#define MIN_OF_THREE(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
#define MAX_OF_THREE(a, b, c) ((a) > (b) ? ((a) > (c) ? (a) : (c)) : ((b) > (c) ? (b) : (c)))
#define AVG_OF_TWO(a, b) (((a)&(b)) + (((a)^(b)) >> 1))

static inline uint64_t hash_64(uint64_t key)
{
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}

#ifndef _PRINT_FORMAT_H_
#define _PRINT_FORMAT_H_

#define NONE "\e[0m" //remove color/font
#define BLACK "\e[0;30m" // black
#define B_BLACK "\e[1;30m" // bold black
#define RED "\e[0;31m" // read
#define B_RED "\e[1;31m" // bold red
#define GREEN "\e[0;32m" // green
#define B_GREEN "\e[1;32m" // bold gren
#define BROWN "\e[0;33m" // brown
#define YELLOW "\e[1;33m" // yellow
#define BLUE "\e[0;34m" // blue
#define B_BLUE "\e[1;34m" // bold blue
#define PURPLE "\e[0;35m" // purple
#define B_PURPLE "\e[1;35m" // bold purple
#define CYAN "\e[0;36m" // cyan
#define B_CYAN "\e[1;36m" // bold cyan
#define GRAY "\e[0;37m" // gray
#define WHITE "\e[1;37m" // white, bold
#define BOLD "\e[1m" // bold
#define UNDERLINE "\e[4m" // underline
#define BLINK "\e[5m" // blink
#define REVERSE "\e[7m" // reverse background and foreground
#define HIDE "\e[8m" // hide
#define STRIKE "\e[9m" // strikethrough
#define CLEAR "\e[2J" // clear
#define CLRLINE "\r\e[K" // clear line

// from https://blog.csdn.net/MoDa_Li/java/article/details/82156888

// #define _RED(string) "\x1b[31m" string "\x1b[0m"
#define _RED(string) "\x1b[91m" string "\x1b[0m"
#define _GREEN(string) "\x1b[32m" string "\x1b[0m"
#define _YELLOW(string) "\x1b[33m" string "\x1b[0m"
#define _BOLD(string) "\x1b[1m" string "\x1b[0m"
#define _BOLD_RED(string) "\x1b[1;91m" string "\x1b[0m"
#define _BOLD_GREEN(string) "\x1b[1;32m" string "\x1b[0m"
#define _BOLD_YELLOW(string) "\x1b[1;33m" string "\x1b[0m"
#define _ERROR_COLOR(string) "\x1b[1;91m" string "\x1b[0m"
#define _WARNING_COLOR(string) "\x1b[1;93m" string "\x1b[0m"
#define _SUCC_COLOR(string) "\x1b[1;32m" string "\x1b[0m"

#define _SUCCESS "\x1b[1;32m"
#define _ERROR "\x1b[1;91m"
#define _WARNING "\x1b[1;93m"

#define WARNING_FORMAT(fp, fmt, ...) fprintf(fp, _WARNING fmt "\x1b[0m", ##__VA_ARGS__)
#define _ERR_WARNING_FORMAT(fmt, ...) fprintf(stderr, _WARNING fmt "\x1b[0m", ##__VA_ARGS__)
#define ERROR_FORMAT(fp, fmt, ...) fprintf(fp, _ERROR fmt "\x1b[0m", ##__VA_ARGS__)
#define _ERR_ERROR_FORMAT(fmt, ...) fprintf(stderr, _ERROR fmt "\x1b[0m", ##__VA_ARGS__)
#define SUCCESS_FORMAT(fp, fmt, ...) fprintf(fp, _SUCCESS fmt "\x1b[0m", ##__VA_ARGS__)
#define _ERR_SUCCESS_FORMAT(fmt, ...) fprintf(stderr, _SUCCESS fmt "\x1b[0m", ##__VA_ARGS__)
#define INFO_FORMAT(fp, fmt, ...) fprintf(fp, fmt, ##__VA_ARGS__)
#define _ERR_INFO_FORMAT(fmt, ...) fprintf(stderr, fmt, ##__VA_ARGS__)

// from https://blog.csdn.net/MoDa_Li/java/article/details/82156888

#endif

#endif
