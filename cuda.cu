#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdlib.h>
#include <stdio.h>

__device__ char mutate_char_cuda(char c)
{
    int mutant = c - 'A';
    mutant = (mutant + 1) % 26;
    return mutant + 'A';
}

__global__ void findBestOffsetKernel(char *dev_main_seq, char *dev_sub_seq, int dev_score_matrix[26][26], int *score_array, int mutant, int main_seq_len, int sub_seq_len)
{
    int offset_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (offset_idx >= main_seq_len - sub_seq_len + 1)
    {
        return;
    }
    int temp_score = 0;

    for (int i = 0; i < sub_seq_len; i++)
    {
        if (i < mutant)
        {
            temp_score += dev_score_matrix[dev_main_seq[offset_idx + i] - 'A'][dev_sub_seq[i] - 'A'];
        }
        else
        {
            temp_score += dev_score_matrix[dev_main_seq[offset_idx + i] - 'A'][mutate_char_cuda(dev_sub_seq[i]) - 'A'];
        }
    }
    score_array[offset_idx] = temp_score;
}

void findBestOffsetCuda(char *dev_main_seq, char *dev_sub_seq, int *dev_score_matrix, int *dev_best_scores, int *best_score, int *best_offset, int *cuda_results, int mutant, int main_seq_len, int sub_seq_len)
{
    *best_score = -2147483648;
    *best_offset = 0;
    int num_offsets = main_seq_len - sub_seq_len + 1;
    int num_blocks = (num_offsets + 256) / 256;
    findBestOffsetKernel<<<num_blocks, 256>>>(dev_main_seq, dev_sub_seq, (int(*)[26])dev_score_matrix, dev_best_scores, mutant, main_seq_len, sub_seq_len);
    cudaDeviceSynchronize();
    cudaMemcpy(cuda_results, dev_best_scores, num_offsets * sizeof(int), cudaMemcpyDeviceToHost);
    for (int i = 0; i < num_offsets; i++)
    {

        if (cuda_results[i] > *best_score)
        {
            *best_score = cuda_results[i];
            *best_offset = i;
        }
    }
}

void cuda_init(char *main_seq, char *sub_seq, int score_matrix[26][26], int main_seq_len, int sub_seq_len, char **dev_main_seq, char **dev_sub_seq, int **dev_score_matrix, int **dev_best_scores)
{
    int num_offsets = main_seq_len - sub_seq_len + 1;
    cudaMalloc((void **)dev_main_seq, main_seq_len * sizeof(char));
    cudaMalloc((void **)dev_sub_seq, sub_seq_len * sizeof(char));
    cudaMalloc((void **)dev_score_matrix, 26 * 26 * sizeof(int));
    cudaMalloc((void **)dev_best_scores, sizeof(int) * num_offsets);

    cudaMemcpy(*dev_main_seq, main_seq, main_seq_len * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(*dev_sub_seq, sub_seq, sub_seq_len * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(*dev_score_matrix, score_matrix, 26 * 26 * sizeof(int), cudaMemcpyHostToDevice);
}

void cuda_free(char *dev_main_seq, char *dev_sub_seq, int *dev_score_matrix, int *dev_best_scores)
{
    cudaFree(dev_main_seq);
    cudaFree(dev_sub_seq);
    cudaFree(dev_score_matrix);
    cudaFree(dev_best_scores);
}