#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "parallel.h"
#include <omp.h>
#include <mpi.h>

int main(int argc, char *argv[])
{

    double cpu_time_used;
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &size); // get number of processes

    parallel(argc, argv, rank, size);

    MPI_Finalize();

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

void findBestOffset(char *main_seq, char *sub_seq, int score_matrix[26][26], int *best_score, int *best_offset, int mutant, int main_seq_len, int sub_seq_len)
{
    int num_offsets = main_seq_len - sub_seq_len + 1;
    *best_score = -2147483648;
    *best_offset = 0;
    int temp_score = 0;
    for (int offset_idx = 0; offset_idx < num_offsets; offset_idx++)
    {
        for (int i = 0; i < sub_seq_len; i++)
        {
            if (i < mutant)
            {
                temp_score += score_matrix[main_seq[offset_idx + i] - 'A'][sub_seq[i] - 'A'];
            }
            else
            {
                temp_score += score_matrix[main_seq[offset_idx + i] - 'A'][mutate_char(sub_seq[i]) - 'A'];
            }
        }
        if (temp_score > *best_score)
        {
            *best_score = temp_score;
            *best_offset = offset_idx;
        }
        temp_score = 0;
    }
}

void findBestScore(char *main_seq, char *sub_seq, int score_matrix[26][26], int *score, int *offset, int *mutant, int rank, int num_procs)
{
    int main_seq_len = strlen(main_seq);
    int sub_seq_len = strlen(sub_seq);
    int num_offsets = main_seq_len - sub_seq_len + 1;
    int num_mutants = sub_seq_len;
    int best_score = -2147483648;
    *score = 0;
    *offset = 0;
    *mutant = 0;
    int mutants_per_proc = num_mutants / num_procs;
    int mutants_start = rank * mutants_per_proc;
    int mutants_end = mutants_start + mutants_per_proc;
    char *dev_main_seq, *dev_sub_seq;
    int *dev_score_matrix, *dev_best_scores, *best_scores;
    cuda_init(main_seq, sub_seq, score_matrix, main_seq_len, sub_seq_len, &dev_main_seq, &dev_sub_seq, &dev_score_matrix, &dev_best_scores);
    best_scores = (int *)malloc(num_offsets * sizeof(int));
    if (rank == num_procs - 1)
    {
        mutants_end = num_mutants;
    }
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int tscore;
        int temp_offset;
        int tbest_offset = 0;
        int tbest_mutant = 0;
        int tbest_score = -2147483648;
#pragma omp for schedule(dynamic, 5)
        for (int mutant_idx = mutants_start; mutant_idx <= mutants_end; mutant_idx++)
        {
            if (tid == 0)
            {
                findBestOffsetCuda(dev_main_seq, dev_sub_seq, dev_score_matrix, dev_best_scores, &tscore, &temp_offset, best_scores, mutant_idx, main_seq_len, sub_seq_len);
            }
            else
            {
                findBestOffset(main_seq, sub_seq, score_matrix, &tscore, &temp_offset, mutant_idx, main_seq_len, sub_seq_len);
            }

            if (tscore > tbest_score)
            {
                tbest_score = tscore;
                tbest_offset = temp_offset;
                tbest_mutant = mutant_idx;
            }
        }
#pragma omp critical
        {
            if (tbest_score > best_score)
            {
                best_score = tbest_score;
                *score = tbest_score;
                *offset = tbest_offset;
                *mutant = tbest_mutant;
            }
        }
    }
    *score = best_score;
    cuda_free(dev_main_seq, dev_sub_seq, dev_score_matrix, dev_best_scores);
    free(best_scores);
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

char mutate_char(char c)
{
    int mutant = c - 'A';
    mutant = (mutant + 1) % 26;
    return mutant + 'A';
}

void mutate_str(char *str, int k)
{
    int len = strlen(str);
    for (int i = 0; i < len; i++)
    {
        if (i < k == 0)
        {
            str[i] = mutate_char(str[i]);
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

void bcastFileAndMatrix(int rank, char *main_seq, char ***sub_sequences, int *sub_sequences_count, int score_matrix[26][26], int argc, char *argv[])
{
    if (rank == 0)
    {
        readInputFile(main_seq, sub_sequences, sub_sequences_count);
        init_score_matrix(score_matrix, argc, argv);
    }
    MPI_Bcast(main_seq, SEQ1_MAX_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(sub_sequences_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(score_matrix, 26 * 26, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0)
    {
        *sub_sequences = (char **)malloc(*sub_sequences_count * sizeof(char *));
        for (int i = 0; i < *sub_sequences_count; i++)
        {
            (*sub_sequences)[i] = (char *)malloc(SEQ2_MAX_LEN * sizeof(char));
        }
    }
    for (int i = 0; i < *sub_sequences_count; i++)
    {
        MPI_Bcast((*sub_sequences)[i], SEQ2_MAX_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
}

void gatherResults(int rank, int num_procs, int *results, int sub_sequences_count)
{
    int *proc_results;
    if (rank == 0)
    {
        proc_results = (int *)malloc(sub_sequences_count * sizeof(int) * 3);
        for (int i = 1; i < num_procs; i++)
        {
            MPI_Recv(proc_results, sub_sequences_count * 3, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < sub_sequences_count; j++)
            {

                if (proc_results[j * 3 + RESULT_SCORE] > results[j * 3 + RESULT_SCORE])
                {
                    results[j * 3 + RESULT_SCORE] = proc_results[j * 3 + RESULT_SCORE];
                    results[j * 3 + RESULT_OFFSET] = proc_results[j * 3 + RESULT_OFFSET];
                    results[j * 3 + RESULT_MUTANT] = proc_results[j * 3 + RESULT_MUTANT];
                }
            }
        }
        free(proc_results);
    }
    else
    {
        MPI_Send(results, sub_sequences_count * 3, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

void parallel(int argc, char *argv[], int rank, int num_procs)
{
    char main_sequence[SEQ1_MAX_LEN];
    char **sub_sequences;
    int score_matrix[26][26];
    int sub_sequences_count;
    bcastFileAndMatrix(rank, main_sequence, &sub_sequences, &sub_sequences_count, score_matrix, argc, argv);
    int *results = (int *)malloc(sub_sequences_count * sizeof(int) * 3);
    double start_time = MPI_Wtime();
    for (int i = 0; i < sub_sequences_count; i++)
    {
        int score;
        int offset;
        int mutant;
        findBestScore(main_sequence, sub_sequences[i], score_matrix, &score, &offset, &mutant, rank, num_procs);
        results[i * 3 + RESULT_SCORE] = score;
        results[i * 3 + RESULT_OFFSET] = offset;
        results[i * 3 + RESULT_MUTANT] = mutant;
    }
    gatherResults(rank, num_procs, results, sub_sequences_count);
    if (rank == 0)
    {
        for (int i = 0; i < sub_sequences_count; i++)
        {
            printf("highest alignment score = %d offset = %d k = %d\n", results[i * 3 + RESULT_SCORE], results[i * 3 + RESULT_OFFSET], results[i * 3 + RESULT_MUTANT]);
        }
    }
    if (rank == 0)
        printf("Time taken: %f\n", MPI_Wtime() - start_time);

    for (int i = 0; i < sub_sequences_count; i++)
    {
        free(sub_sequences[i]);
    }
    free(sub_sequences);
}