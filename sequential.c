#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "sequential.h"

int main(int argc, char *argv[])
{
    // start time
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    sequential(argc, argv);

    // end time
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken: %f\n", cpu_time_used);

    return 0;
}

void readInputFile(char *main_sequence, char ***sub_sequences, int *sub_sequences_count)
{
    FILE *file = stdin;
    if (file == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    fscanf(file, "%s", main_sequence);
    strToUpper(main_sequence);
    fscanf(file, "%d", sub_sequences_count);
    *sub_sequences = (char **)malloc(*sub_sequences_count * sizeof(char *));
    for (int i = 0; i < *sub_sequences_count; i++)
    {
        (*sub_sequences)[i] = (char *)malloc(SEQ2_MAX_LEN * sizeof(char));
        fscanf(file, "%s", (*sub_sequences)[i]);
        strToUpper((*sub_sequences)[i]);
    }
    fclose(file);
}

void findBestScore(char *main_seq, char *sub_seq, int score_matrix[26][26], int *score, int *offset, int *mutant)
{
    int main_seq_len = strlen(main_seq);
    int sub_seq_len = strlen(sub_seq);
    int num_offsets = main_seq_len - sub_seq_len + 1;
    int num_mutants = sub_seq_len;
    int best_score = -2147483648;
    *score = 0;
    *offset = 0;
    *mutant = 0;
    for (int mutant_idx = 0; mutant_idx <= num_mutants; mutant_idx++)
    {
        for (int offset_idx = 0; offset_idx < num_offsets; offset_idx++)
        {
            for (int i = 0; i < sub_seq_len; i++)
            {
                if (i < mutant_idx)
                {
                    *score += score_matrix[main_seq[offset_idx + i] - 'A'][sub_seq[i] - 'A'];
                }
                else
                {
                    *score += score_matrix[main_seq[offset_idx + i] - 'A'][mutant_char(sub_seq[i]) - 'A'];
                }
            }
            if (*score > best_score)
            {
                best_score = *score;
                *offset = offset_idx;
                *mutant = mutant_idx;
            }
            *score = 0;
        }
    }
    *score = best_score;
}

char toUpper(char c)
{
    if (c >= 'a' && c <= 'z')
    {
        c = c - 'a' + 'A';
    }
    return c;
}

void strToUpper(char *str)
{
    int len = strlen(str);
    for (int i = 0; i < len; i++)
    {
        str[i] = toUpper(str[i]);
    }
}

char mutant_char(char c)
{
    int mutant = c - 'A';
    mutant = (mutant + 1) % 26;
    return mutant + 'A';
}

void mutant_str(char *str, int k)
{
    int len = strlen(str);
    for (int i = 0; i < len; i++)
    {
        if (i < k == 0)
        {
            str[i] = mutant_char(str[i]);
        }
    }
}

void init_score_matrix(int score_matrix[26][26], int argc, char *argv[])
{
    memset(score_matrix, 0, sizeof(int) * 26 * 26);
    if (argc > 1)
    {
        FILE *file = fopen(argv[1], "r");
        if (file == NULL)
        {
            printf("Error opening file!\n");
            exit(1);
        }
        for (int i = 0; i < 26; i++)
        {
            for (int j = 0; j < 26; j++)
            {
                fscanf(file, "%d", &score_matrix[i][j]);
            }
        }
        fclose(file);
    }
    else
    {
        for (int i = 0; i < 26; i++)
        {
            score_matrix[i][i] = 1;
        }
    }
}

void sequential(int argc, char *argv[])
{
    char main_sequence[SEQ1_MAX_LEN];
    char **sub_sequences;
    int score_matrix[26][26];
    int sub_sequences_count;
    readInputFile(main_sequence, &sub_sequences, &sub_sequences_count);
    init_score_matrix(score_matrix, argc, argv);
    int num_mutants = 0;
    int total_score = 0;
    for (int i = 0; i < sub_sequences_count; i++)
    {
        int score;
        int offset;
        int mutant;
        findBestScore(main_sequence, sub_sequences[i], score_matrix, &score, &offset, &mutant);
        printf("highest alignment score = %d offset = %d k = %d\n", score, offset, mutant);
        free(sub_sequences[i]);
    }
    free(sub_sequences);
}