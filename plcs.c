/*based on Jiang Yanyan's course miniLab, change a little bit 
usage: ./plcs < input.txt, output should be 2820*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "thread.h"
#include "thread-sync.h"

#define MAXN 10000
#define num_threads 16
#define num_block 100
#define block_size 100

int T, N, M;
char A[MAXN + 1];
char B[MAXN + 1];
int dp[MAXN][MAXN];
int result;
spinlock_t *lk;
#define DP(x, y) (((x) >= 0 && (y) >= 0) ? dp[x][y] : 0)
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MAX3(x, y, z) MAX(MAX(x, y), z)

int done[num_block][num_block] = {0};
int stack[4 * num_block + 4];
int stack_top = 0;
void stack_push(int m, int n)
{
    stack[stack_top] = m;
    stack[stack_top + 1] = n;
    stack_top += 2;

    // printf("push(%d,%d)\n", m, n);
}
void stack_pop(int *m, int *n)
{
    *m = stack[stack_top - 2];
    *n = stack[stack_top - 1];
    stack_top -= 2;

    // printf("pop(%d,%d)\n", *m, *n);
}
void fill_block(int m, int n)
{
    for (int i = m * block_size; i < (m + 1) * block_size; i++)
    {
        for (int j = n * block_size; j < (n + 1) * block_size; j++)
        {
            int skip_a = DP(i - 1, j);
            int skip_b = DP(i, j - 1);
            int take_both = DP(i - 1, j - 1) + (A[i] == B[j]);
            dp[i][j] = MAX3(skip_a, skip_b, take_both);
        }
    }
    // printf("done(%d,%d)\n", m, n);
}
void Tworker(int id)
{
    int m, n;
    while (1)
    {
        spin_lock(lk);
        if (done[num_block - 1][num_block - 1] == 1)
        {
            spin_unlock(lk);
            return;
        }
        else if (stack_top != 0)
        {
            stack_pop(&m, &n);
            spin_unlock(lk);
            fill_block(m, n);
            spin_lock(lk);
            done[m][n] = 1;
            if (m == 0 && n < num_block - 1)
            {
                stack_push(m, n + 1);
            }
            if (n == 0 && m < num_block - 1)
            {
                stack_push(m + 1, n);
            }
            if (m > 0 && n < num_block - 1 && done[m - 1][n + 1] == 1)
            {
                stack_push(m, n + 1);
            }
            if (n > 0 && m < num_block - 1 && done[m + 1][n - 1] == 1)
            {
                stack_push(m + 1, n);
            }
            spin_unlock(lk);
        }
        else
        {
            spin_unlock(lk);
        }
    }
}

int main(int argc, char *argv[])
{
    assert(scanf("%s%s", A, B) == 2);
    lk = (int *)malloc(sizeof(int));
    *lk = SPIN_INIT();
    // Add preprocessing code here
    stack_push(0, 0);
    for (int i = 0; i < num_threads; i++)
    {
        create(Tworker);
    }
    join(); // Wait for all workers
    result = dp[MAXN - 1][MAXN - 1];
    printf("%d\n", result);
}
