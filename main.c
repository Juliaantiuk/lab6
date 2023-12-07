//THIS IS FINAL
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <conio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#define MAX_NUM 15
#define MAX_PRECIS 15
#define MAX_COEF 1e+6
#define MIN_COEF 1e-6

#define RED "\x1b[31m"
#define BLUE "\x1b[36m"
#define RESET "\x1b[0m"
#define GREEN "\033[0;32m"
#define ESC 27

void prnt_greeting() {
    printf(BLUE"Hello! This is a program for solving the system of linear algebraic equations with a simple iteration method\n"RESET);
    printf(BLUE"The maximum number of equations is %d, the maximum precision is %d, the maximum coeficient is |%.e|, the minimum is |%.e| or 0\n"RESET, MAX_NUM, MAX_PRECIS, MAX_COEF, MIN_COEF);
    printf(BLUE"An example of SLAE:\n"RESET);
    printf(BLUE"a[1][1]*x[1] + a[1][2]*x[2] + ... + a[1][n]*x[n] = b[1]\n"RESET);
    printf(BLUE"a[2][1]*x[1] + a[2][2]*x[2] + ... + a[2][n]*x[n] = b[2]\n...\n"RESET);
    printf(BLUE"a[n][1]*x[1] + a[n][2]*x[2] + ... + a[n][n]*x[n] = b[n]\n\n"RESET);
}

bool is_esc() {
    printf(BLUE"Press ESC key to exit, or any other key to try again: "RESET);
    char escChoice = _getch();
    printf("\n");
    system("cls");
    return escChoice == ESC;
}

bool is_input_valid(void* inp_ptr, const char* format) {
    int start_space = 0;
    bool isInputValid = false;
    int newline = 0;
    int result = scanf(format, &start_space, inp_ptr, &newline);
    rewind(stdin);
    if (start_space == 0 && result == 2 && newline == '\n') {
        isInputValid = true;
    }
    else {
        printf(RED"Error! Invalid input format!\n\n"RESET);
    }
    return isInputValid;
}

bool is_restriction_valid(int num, int low_restr, int up_restr) {
    if (!(num >= low_restr && num <= up_restr)) {
        printf(RED"Error! The number must be between %d and %d\n"RESET, low_restr, up_restr);
        return false;
    }
    return true;
}

bool is_coef_valid(double num, double low_rest, double up_restr) {
    if (!((fabs(num) >= low_rest && fabs(num) <= up_restr) || num == 0)) {
        printf(RED"Error! The module of coeficient must be between %.e and %.e or 0\n"RESET, low_rest, up_restr);
        return false;
    }
    return true;
}

bool is_equation_valid(double* a_line, int n, int index) {
    double sum = 0;
    double main = a_line[index];
    for (int i = 0; i < n; i++) {
        if (i != index) {
            sum += fabs(a_line[i]);
        }
    }
    if (fabs(main) < sum) {
        printf(RED"\nThe system with such coeficients is not convergent! Try again!\n\n"RESET);
        return false;
    }
    return true;
}

void take_inp(int num, double** main_m, double* free_m) {
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < num; j++) {
            do {
                printf(BLUE"Enter a[%d][%d]: "RESET, i + 1, j + 1);
            } while (!is_input_valid(&main_m[i][j], " %n%lf%c") || !is_coef_valid(main_m[i][j], MIN_COEF, MAX_COEF));
        }
        if (!is_equation_valid(main_m[i], num, i)) {
            i--;
        }
        printf("\n");
    }

    for (int i = 0; i < num; i++) {
        do {
            printf(BLUE"Enter b[%d]: "RESET, i + 1);
        } while (!is_input_valid(&free_m[i], " %n%lf%c") || !is_coef_valid(free_m[i], MIN_COEF, MAX_COEF));
    }
    printf("\n");
}

double convert_eps(int eps) {
    double correction = (eps - 1) / pow(10, eps);
    double epsilon = (eps / pow(10, eps)) - correction;
    return epsilon;
}

void allocate_memory(double*** a, double** b, double** x, double** xp, double** delta, int n) {
    *a = (double**)calloc(n, sizeof(double*));
    *b = (double*)calloc(n, sizeof(double));
    *x = (double*)calloc(n, sizeof(double));
    *xp = (double*)calloc(n, sizeof(double));
    *delta = (double*)calloc(n, sizeof(double));

    for (int i = 0; i < n; i++) {
        (*a)[i] = (double*)calloc(n, sizeof(double));
    }
}


void free_memory(double** a, double* b, double* x, double* xp, int n) {
    for (int i = 0; i < n; i++) {
        free(a[i]);
    }
    free(a);
    free(b);
    free(x);
    free(xp);
}

bool is_zero(double** main_matr, double* free_matr, int n) {
    bool isAllMainZero = true;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (main_matr[i][j] != 0) {
                isAllMainZero = false;
            }
        }
        if (isAllMainZero) {
            printf(GREEN"This system does not have roots or has an infinite amount of roots\n"RESET);
            return true;
        }
    }
    return false;
}

void solve_equation(double** main_matr, double* free_matrix, double* x, double* x_prev, double* delta, int n, double eps) {
    double epsilon = convert_eps(eps);
    double max_delta = 0, sum = 0;

    do {
        max_delta = 0;
        for (int i = 0; i < n; i++) {
            x_prev[i] = x[i];
            sum = 0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum += (main_matr[i][j] * x[j]);
                }
            }
            x[i] = (free_matrix[i] - sum) / main_matr[i][i];
            delta[i] = fabs((x[i] - x_prev[i]) / x[i]);
            if (delta[i] > max_delta) {
                max_delta = delta[i];
            }
        }
    } while (max_delta >= epsilon);
}

void check(int num, double** main_matr, double* x_matr, int digits) {
    double res = 0;
    for (int i = 0; i < num; i++) {
        printf(BLUE"%lf * %.*lf"RESET, main_matr[i][0], digits, x_matr[0]);
        res += main_matr[i][0] * x_matr[0];
        for (int j = 1; j < num; j++) {
            printf(BLUE" + %lf * %.*lf"RESET, main_matr[i][j], digits, x_matr[j]);
            res += main_matr[i][j] * x_matr[j];
        }
        printf(GREEN" = %.*lf\n"RESET, digits, res);
        res = 0;
    }
}

int main() {
    int n = 0, e = 0;
    do {
        prnt_greeting();
        do {
            printf(BLUE"Enter the number of equations (2 - %d): "RESET, MAX_NUM);
        } while (!is_input_valid(&n, " %n%d%c") || !is_restriction_valid(n, 2, MAX_NUM));
        do {
            printf(BLUE"Enter the presicion (max: %d): "RESET, MAX_PRECIS);
        } while (!is_input_valid(&e, " %n%d%c") || !is_restriction_valid(e, 1, MAX_PRECIS));
        printf("\n");
        double** matr_a = 0, * matr_b = 0, * x = 0, * prev_x = 0, * delta = 0;
        allocate_memory(&matr_a, &matr_b, &x, &prev_x, &delta, n);
        take_inp(n, matr_a, matr_b);
        if (!is_zero(matr_a, matr_b, n)) {
            solve_equation(matr_a, matr_b, x, prev_x, delta, n, e);
            for (int i = 0; i < n; i++) {
                printf("%sx[%d] = %s%s%.*lf%s\n", BLUE, i + 1, RESET, GREEN, e, x[i], RESET);
            }
            printf("\n");
            check(n, matr_a, x, e);
        }
        free_memory(matr_a, matr_b, x, prev_x, n);
        printf("\n");
    } while (!is_esc());
    return 0;
}
