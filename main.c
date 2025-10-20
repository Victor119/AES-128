#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "sboxMatrix.h"

#define _GNU_SOURCE 
#define Nb 4 // The number of columns for a state in AES
#define Nr 10 // Number of rounds for AES-128
#define max(a,b) ((a) > (b) ? (a) : (b))


typedef struct Node {
    int power;
    int coeff;
    struct Node* next;
} Node;

typedef struct {
    int first;  // exponent
    int second; // coeficient
} Pair;

char* binar_la_polinom(const char* binar) {
    int len = strlen(binar);
    if (len == 0) return NULL;
    
    char* rezultat = (char*)malloc(1000 * sizeof(char));
    rezultat[0] = '\0';
    
    int primul = 1;
    
    for (int i = 0; i < len; i++) {
        if (binar[i] == '1') {
            int grad = len - 1 - i;
            
            if (primul == 0) {
                strcat(rezultat, " + ");
            }
            primul = 0;
            
            if (grad == 0) {
                strcat(rezultat, "1");
            } else if (grad == 1) {
                strcat(rezultat, "x");
            } else {
                char termen[50];
                sprintf(termen, "x^%d", grad);
                strcat(rezultat, termen);
            }
        }
    }
    
    if (primul) {
        strcpy(rezultat, "0");
    }
    
    return rezultat;
}

unsigned char skip_spaces(char *p, char **pp) {
    unsigned char c;
    while (isspace(c = *p))
        p++;
    *pp = p;
    return c;
}


void show_error(const char *error, const char *str, const char *p) {
    fprintf(stderr, "error: %s\n%s%*s^\n", error, str, (int)(p - str), "");
}

Pair* parse_polinom(const char* polinom, int* num_terms) {
    #define TERM_COUNT 100
    static Pair perechi[TERM_COUNT];
    int coef[TERM_COUNT] = { 0 };
    int expo[TERM_COUNT] = { 0 };
    int counter;
    char *p = (char*)polinom;

    for (counter = 0; counter < TERM_COUNT; counter++) {
        int sign = 1;
        unsigned char c = skip_spaces(p, &p);
        if (c == '\0')
            break;
        if (c == '+') {
            c = skip_spaces(p + 1, &p);
        } else if (c == '-') {
            sign = -1;
            c = skip_spaces(p + 1, &p);
        }
        if (c == '\0') {
            show_error("missing term", polinom, p);
            break;
        }
        if (isdigit(c)) {
            coef[counter] = sign * strtol(p, &p, 10);
            c = skip_spaces(p, &p);
        } else {
            coef[counter] = sign;
        }
        if (c == 'x') {
            c = skip_spaces(p + 1, &p);
            if (c == '^') {
                c = skip_spaces(p + 1, &p);
                if (isdigit(c)) {
                    expo[counter] = strtol(p, &p, 10);
                    c = skip_spaces(p, &p);
                } else {
                    show_error("missing exponent", polinom, p);
                    break;
                }
            } else {
                expo[counter] = 1;
            }
        } else {
            expo[counter] = 0;
        }
        if (c != '\0' && c != '+' && c != '-') {
            show_error("syntax error", polinom, p);
            break;
        }
    }

    // Creăm vectorul de perechi
    for (int i = 0; i < counter; i++) {
        perechi[i].first = expo[i];
        perechi[i].second = coef[i];
    }
    
    *num_terms = counter;
    return perechi;
}

Node* creeazaListaDinPerechi(Pair* perechi, int num_terms) {
    if (num_terms == 0) return NULL;

    // primul nod
    Node* head = (Node*)malloc(sizeof(Node));
    head->power = perechi[0].first;
    head->coeff = perechi[0].second;
    head->next = NULL;

    Node* current = head;

    // adăugăm restul nodurilor
    for (int i = 1; i < num_terms; i++) {
        Node* nou = (Node*)malloc(sizeof(Node));
        nou->power = perechi[i].first;
        nou->coeff = perechi[i].second;
        nou->next = NULL;

        current->next = nou;
        current = nou;
    }

    return head;
}

void afiseazaLista(Node* head) {
    Node* temp = head;
    while (temp != NULL) {
        printf("(%d, %d) ", temp->power, temp->coeff);
        temp = temp->next;
    }
    printf("\n");
}

void stergeLista(Node* head) {
    Node* temp;
    while (head != NULL) {
        temp = head;
        head = head->next;
        free(temp);
    }
}

Node *multiplyPolynomials(Node *poly1, Node *poly2) {
    if (poly1 == NULL || poly2 == NULL) return NULL;

    // Gasim puterea maxima pentru a aloca array-ul
    int maxPower = 0;
    for (Node *curr1 = poly1; curr1 != NULL; curr1 = curr1->next) {
        for (Node *curr2 = poly2; curr2 != NULL; curr2 = curr2->next) {
            int power = curr1->power + curr2->power;
            if (power > maxPower) maxPower = power;
        }
    }

    // Alocam array pentru coeficienti
    int *coeffs = (int*)calloc(maxPower + 1, sizeof(int));
    if (coeffs == NULL) return NULL;

    // Multiply each node of poly2 with poly1
    for (Node *curr1 = poly1; curr1 != NULL; curr1 = curr1->next) {
        for (Node *curr2 = poly2; curr2 != NULL; curr2 = curr2->next) {
            int power = curr1->power + curr2->power;
            int coeff = curr1->coeff * curr2->coeff;
            coeffs[power] += coeff;
        }
    }

    // Create the result linked list
    Node *resultHead = NULL;
    Node *prev = NULL;

    for (int it = maxPower; it >= 0; it--) {
        if (coeffs[it] != 0) { // Adaugam doar termenii nenuli
            Node *newNode = (Node*)malloc(sizeof(Node));
            newNode->power = it;
            newNode->coeff = coeffs[it];
            newNode->next = NULL;

            if (resultHead == NULL) {
                resultHead = newNode;
                prev = newNode;
            }
            else {
                prev->next = newNode;
                prev = newNode;
            }
        }
    }

    free(coeffs);
    return resultHead;
}

void printPolynomial(Node *head) {
    while (head != NULL) {
        printf("%dx^%d", head->coeff, head->power);
        if (head->next != NULL)
            printf(" + ");
        head = head->next;
    }
    printf("\n");
}

char* polynomialToString(Node* head) {
    if (!head) return strdup("0");

    char buffer[2048];
    buffer[0] = '\0';
    Node* temp = head;
    int first = 1;

    while (temp != NULL) {
        char term[64];
        int coeff = temp->coeff;
        int power = temp->power;

        if (coeff == 0) {
            temp = temp->next;
            continue;
        }

        if (!first && coeff > 0)
            strcat(buffer, " + ");
        else if (coeff < 0)
            strcat(buffer, " - ");

        if (abs(coeff) != 1 || power == 0)
            sprintf(term, "%d", abs(coeff));
        else
            term[0] = '\0';

        if (power > 0) {
            if (term[0] != '\0')
                strcat(term, "x");
            else
                strcpy(term, "x");
        }
        if (power > 1) {
            char exp_str[16];
            sprintf(exp_str, "^%d", power);
            strcat(term, exp_str);
        }

        strcat(buffer, term);
        first = 0;
        temp = temp->next;
    }

    return strdup(buffer);
}

char* inlocuieste_x8(const char* polinom) {
    // Inlocuieste orice "x^8" sau "1x^8" cu "x^4 + x^3 + x + 1"
    const char* pattern1 = "1x^8";
    const char* pattern2 = "x^8";
    const char* replacement = "x^4 + x^3 + x + 1";

    char* result = malloc(4096);
    result[0] = '\0';

    const char* p = polinom;
    while (*p) {
        if (strncmp(p, pattern1, strlen(pattern1)) == 0) {
            strcat(result, replacement);
            p += strlen(pattern1);
        } else if (strncmp(p, pattern2, strlen(pattern2)) == 0) {
            strcat(result, replacement);
            p += strlen(pattern2);
        } else {
            strncat(result, p, 1);
            p++;
        }
    }

    return result;
}

char* polynomialPairsToBinary(Pair* pairs, int num_terms) {
    // gasim exponentul maxim pentru a sti lungimea binarului
    int max_power = pairs[0].first;  // primul are exponentul maxim

    char* binary = (char*)malloc(9); // 8 biti + terminator
    if (!binary) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    // initializam toti bitii cu 0
    for (int i = 0; i < 8; i++)
        binary[i] = '0';
    binary[8] = '\0';

    int counter = 0;
    // setam bitii corespunzători termenilor existenti
    for (int i = 0; i < num_terms; i++) {
        int exponent = pairs[i].first;
        if (exponent >= 0 && exponent < 8) {
            int position = 7 - exponent;
            binary[position] = '1';
        }
    }

    return binary;
}

char* multiplication_in_GF_2_power_8(const char* userInput, const char* userInput2) {
    int num_terms;
    Pair* perechi;

    char* polinom1 = binar_la_polinom(userInput);

    if (polinom1) {
        polinom1[strcspn(polinom1, "\n")] = 0;
        perechi = parse_polinom(polinom1, &num_terms);

        printf("Vectorul de perechi (exponent, coeficient): %s \n", polinom1);
        for (int i = 0; i < num_terms; i++) {
            printf("(%d, %d) ", perechi[i].first, perechi[i].second);
        }
        printf("\n");

        free(polinom1);
    }
    Node* poly1 = creeazaListaDinPerechi(perechi, num_terms);

    char* polinom2 = binar_la_polinom(userInput2);

    if (polinom2) {
        polinom2[strcspn(polinom2, "\n")] = 0;
        perechi = parse_polinom(polinom2, &num_terms);

        printf("Vectorul de perechi (exponent, coeficient): %s \n", polinom2);
        for (int i = 0; i < num_terms; i++) {
            printf("(%d, %d) ", perechi[i].first, perechi[i].second);
        }
        printf("\n");

        free(polinom2);
    }
    Node* poly2 = creeazaListaDinPerechi(perechi, num_terms);

    Node *result = multiplyPolynomials(poly1, poly2);
    printPolynomial(result);

    // curatam memoria
    stergeLista(poly1);
    stergeLista(poly2);

    char* poly_str = polynomialToString(result);
    printf("\nPolinom rezultat (ca sir): %s\n", poly_str);

    char* binar_final = NULL;

    if (strstr(poly_str, "x^8") != NULL) {
        char* reduced_poly = inlocuieste_x8(poly_str);
        printf("Dupa inlocuirea lui x^8: %s\n", reduced_poly);

        // textul pe care vrem să-l extragem
        const char* pattern = "x^4 + x^3 + x + 1 + ";
        size_t pattern_len = strlen(pattern);

        if (strncmp(reduced_poly, pattern, pattern_len) == 0) {
            // Extragem partea ce ramane după eliminarea pattern-ului
            char* remaining = reduced_poly + pattern_len;

            printf("Partea extrasa (eliminata): %s\n", pattern);
            printf("Dupa scoaterea lui x^8: %s\n", remaining);

            char* final_poly = malloc(strlen(remaining) + 10);
            strcpy(final_poly, remaining);
            strcat(final_poly," + 1");
            printf("Polinom final: %s\n", final_poly);
            
            Pair* final_poly_pairs = parse_polinom(final_poly, &num_terms);

            printf("Vectorul de perechi (exponent, coeficient):\n");
            for (int i = 0; i < num_terms; i++) {
                printf("(%d, %d) ", final_poly_pairs[i].first, final_poly_pairs[i].second);
            }
            printf("\n");

            int new_num_terms = 0;
            for (int i = 0; i < num_terms; i++) {
                // Daca coeficientul este impar, il pastram
                if (final_poly_pairs[i].second % 2 != 0) {
                    final_poly_pairs[new_num_terms].first = final_poly_pairs[i].first;

                    // Daca coeficientul este impar si mai mare ca 1, il facem 1
                    if (final_poly_pairs[i].second > 1) {
                        final_poly_pairs[new_num_terms].second = 1;
                    } 
                    // Daca coeficientul este impar si mai mic ca -1, il facem -1
                    else if (final_poly_pairs[i].second < -1) {
                        final_poly_pairs[new_num_terms].second = -1;
                    }
                    // Altfel (coeficientul este 1 sau -1), il pastram asa
                    else {
                        final_poly_pairs[new_num_terms].second = final_poly_pairs[i].second;
                    }
                    
                    new_num_terms++;
                }
                // Daca coeficientul este par, nu il mai adaugam (il eliminam)
            }

            // Actualizam numarul de termeni
            num_terms = new_num_terms;

            printf("Dupa procesare (eliminare coeficienti pari si normalizare impari):\n"); 
            for (int i = 0; i < num_terms; i++) {
                printf("(%d, %d) ", final_poly_pairs[i].first, final_poly_pairs[i].second);
            }
            printf("\n");

            binar_final = polynomialPairsToBinary(final_poly_pairs, num_terms);
            printf("Forma binara (8 biti) a polinomului final: %s\n", binar_final);

            free(final_poly);
            free(reduced_poly);
        }
    }
    else{
        int new_num_terms = 0;
        Pair* final_poly_pairs = parse_polinom(poly_str, &num_terms);

        for (int i = 0; i < num_terms; i++) {
            // Daca coeficientul este impar, il pastram
            if (final_poly_pairs[i].second % 2 != 0) {
                final_poly_pairs[new_num_terms].first = final_poly_pairs[i].first;

                // Daca coeficientul este impar si mai mare ca 1, il facem 1
                if (final_poly_pairs[i].second > 1) {
                    final_poly_pairs[new_num_terms].second = 1;
                } 
                // Daca coeficientul este impar si mai mic ca -1, il facem -1
                else if (final_poly_pairs[i].second < -1) {
                    final_poly_pairs[new_num_terms].second = -1;
                }
                // Altfel (coeficientul este 1 sau -1), il pastram asa
                else {
                    final_poly_pairs[new_num_terms].second = final_poly_pairs[i].second;
                }
                    
                new_num_terms++;
            }
            // Daca coeficientul este par, nu il mai adaugam (il eliminam)
        }

        // Actualizam numarul de termeni
        num_terms = new_num_terms;

        printf("Dupa procesare (eliminare coeficienti pari si normalizare impari):\n"); 
        for (int i = 0; i < num_terms; i++) {
            printf("(%d, %d) ", final_poly_pairs[i].first, final_poly_pairs[i].second);
        }
        printf("\n");

        binar_final = polynomialPairsToBinary(final_poly_pairs, num_terms);
        printf("Forma binara (8 biti) a polinomului final: %s\n", binar_final);
    }

    free(poly_str);
    stergeLista(result);

    return binar_final;
}

unsigned createMask(unsigned a, unsigned b)
{
    unsigned r = 0;
    for (unsigned i=a; i<=b; i++)
        r |= 1 << i;

    return r;
}

int BinaryStringToInt(const char* binaryCharArray) {
    int total = 0;
    for (size_t i = 0; i < strlen(binaryCharArray); i++) {
        total = (total << 1) | (binaryCharArray[i] - '0');
    }
    return total;
}

// function to reverse the string
void reverse(char *bin, int left, int right) {
    while (left < right) {
        char temp = bin[left];
        bin[left] = bin[right];
        bin[right] = temp;
        left++;
        right--;
    }
}

// function to convert decimal to binary
char* decToBinary(int n) {
    char* bin = (char*) malloc(9); // 8 biți + terminator
    bin[8] = '\0';
    for(int i = 7; i >= 0; i--) {
        bin[i] = (n & 1) ? '1' : '0';
        n >>= 1;
    }
    return bin;
}

void xor_two_matrix_columns(char*** matrix, int col1, int col2, char** result){
    for (int i = 0; i < 4; i++) {
        if (!matrix[i][col1] || !matrix[i][col2]) {
            result[i] = decToBinary(0);
        } else {
            int val1 = BinaryStringToInt(matrix[i][col1]);
            int val2 = BinaryStringToInt(matrix[i][col2]);
            int xor_result = val1 ^ val2;
            result[i] = decToBinary(xor_result);
        }
    }
}


char* binary_string(const char *hex_string) {
    if (!hex_string) return NULL;

    size_t len = strlen(hex_string);
    char *binary = (char*)malloc(len * 4 + 1); // 4 biți per hex + terminator '\0'
    if (!binary) return NULL;

    binary[0] = '\0'; // inițializează șirul gol

    for (size_t i = 0; i < len; i++) {
        char c = hex_string[i];
        switch (c) {
            case '0': strcat(binary, "0000"); break;
            case '1': strcat(binary, "0001"); break;
            case '2': strcat(binary, "0010"); break;
            case '3': strcat(binary, "0011"); break;
            case '4': strcat(binary, "0100"); break;
            case '5': strcat(binary, "0101"); break;
            case '6': strcat(binary, "0110"); break;
            case '7': strcat(binary, "0111"); break;
            case '8': strcat(binary, "1000"); break;
            case '9': strcat(binary, "1001"); break;
            case 'A': case 'a': strcat(binary, "1010"); break;
            case 'B': case 'b': strcat(binary, "1011"); break;
            case 'C': case 'c': strcat(binary, "1100"); break;
            case 'D': case 'd': strcat(binary, "1101"); break;
            case 'E': case 'e': strcat(binary, "1110"); break;
            case 'F': case 'f': strcat(binary, "1111"); break;
            case '.': strcat(binary, "."); break;
            default:
                free(binary);
                return NULL; // caracter invalid
        }
    }

    return binary;
}

// Function to add two Binary Strings
char* sum(char *b1, char *b2, int l1, int l2)
{
    int carry = 0, temp, num1, num2, i;

    char *result = (char*) malloc( max(l1, l2) + 2 );
    if (!result) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    result[max(l1, l2) + 1] = '\0';

    // This loop will add Both Strings till
    // second have characters in it.
    while (l2 > 0) {
        num1 = b1[l1 - 1] - '0';
        num2 = b2[l2 - 1] - '0';
        temp = num1 + num2 + carry;
        carry = 0;
        if (temp >= 2) {
            carry = 1;
            temp = temp % 2;
        }
        result[l1] = temp + '0';
        l1--;
        l2--;
    }

    // This loop will add leftover
    // characters of first Strings.
    while (l1 - 1 >= 0) {
        temp = b1[l1 - 1] + carry - '0';
        if (temp >= 2) {
            carry = 1;
            temp = temp % 2;
        }
        result[l1] = temp + '0';
        l1--;
    }

    // Add last carry to result string
    if (carry) {
        result[0] = '1';
    }
    else {
        // if there is no carry then we will shift
        // each character by one place in left side.
        for (i = 0; i < strlen(result) - 1; i++) {
            result[i] = result[i + 1];
        }
        result[strlen(result) - 1] = '\0';
    }

    return result;
}

char* xoring(char* a, char* b, int n) {
    char* ans = malloc(n + 1);
    for (int i = 0; i < n; i++)
        ans[i] = (a[i] == b[i]) ? '0' : '1';
    ans[n] = '\0';
    return ans;
}




void AddRoundKey(char*** binary_state, char*** w2, uint8_t sbox_matrix[16][16], int expanded_cols, int rows) {
    int i, j, k;
    
    printf("\n");
    for(j = 4; j < expanded_cols; j++) {
        char* temp_col[4] = {NULL, NULL, NULL, NULL};

        if((j+1)%4 != 0){
            xor_two_matrix_columns(w2, j-4, j-1, temp_col);

            for(i = 0; i < rows; i++) {
                if (temp_col[i]) {
                    w2[i][j] = temp_col[i];
                    printf("%s ", w2[i][j]);
                }
            }
            printf("\n");
        }
        if((j+1)%4 == 0 && j>4){
            char* T_mat[4] = {w2[1][j-1], w2[2][j-1], w2[3][j-1], w2[0][j-1]};
            for (k = 0; k < 4; k++){
                char *source = T_mat[k];
                if (!source) source = "00000000";

                char *buf1 = malloc(5);
                char *buf2 = malloc(5);
                if (!buf1 || !buf2) {
                    perror("malloc");
                    exit(EXIT_FAILURE);
                }
                buf1[0] = '\0';
                buf2[0] = '\0';

                //char tmp8[9]; tmp8[8]='\0';
                //strncpy(tmp8, source, 8);
                //tmp8[8] = '\0';

                strncat(buf1, source, 4);      // primele 4 biți
                strncat(buf2, source + 4, 4);  // următorii 4 biți

                int row2 = BinaryStringToInt(buf1);
                int col2 = BinaryStringToInt(buf2);

                free(buf1);
                free(buf2);
                    
                uint8_t val = sbox_matrix[row2][col2];
                w2[k][j] = decToBinary(val);

                
            }
            for (i = 0; i < 4; ++i) {
                if (w2[i][j]) printf("%s ", w2[i][j]);
            }
            printf("\n");

            //de convertit r_i to string in binary form
            int factor = ((j+1-4)/4);

            char *x = "00000010";
            char *r_i = "00000010";

            size_t size_x = strlen(x);
            size_t size_r_i = strlen(r_i);

            for(i=0; i<factor; i++){
                r_i = sum(r_i, x, size_x, size_r_i);
                size_r_i = strlen(r_i);
            }

            // ne asiguram ca r_i are exact 8 biti
            while(strlen(r_i) < 8) {
                char *temp = malloc(9);
                sprintf(temp, "0%s", r_i);
                free(r_i);
                r_i = temp;
            }

            int bin_val = BinaryStringToInt(w2[0][j]);
            char *bin_char = decToBinary(bin_val);

            printf("e=%s ", w2[0][j]);
            printf("\n");
            printf("r_i=%s ", r_i);
            printf("\n");
            
            w2[0][j] = xoring(bin_char, r_i, 8);
            free(bin_char);

            
            printf("e xor r_i=%s ", w2[0][j]);
            printf("\n");
            printf("T(e xor r_i, f, g, h) = ");
            for (i = 0; i < 4; ++i) {
                if (w2[i][j]) printf("%s ", w2[i][j]);
            }
            
            printf("\n");

            printf("w_8 = x_4 xor T(e xor r_i, f, g, h) :");
            //T(e xor r_i, f, g, h) se afla pe coloana 7 din w2 
            xor_two_matrix_columns(w2, j-4, j, temp_col);

            for(i = 0; i < rows; i++) {
                if (temp_col[i]) {
                    w2[i][j] = temp_col[i];
                    printf("%s ", w2[i][j]);
                }
            }
            printf("\n");
        }
    }

    //binary_state xor w2[4..7]
    for (i = 0; i < 4; i++) {
        for(j=0; j<4; j++){
            int val1 = BinaryStringToInt(binary_state[j][i]);
            int val2 = BinaryStringToInt(w2[j][i+4]);
            int xor_result = val1 ^ val2;
            binary_state[j][i] = decToBinary(xor_result);
        }
    }
}

void SubBytes(char*** binary_state, uint8_t sbox_matrix[16][16]) {
    int j, k;
    
    printf("After SubBytes:\n");
    
    for(j = 0; j < 4; j++) {
        for(k = 0; k < 4; k++) {
            char *source = binary_state[j][k];
            
            char *buf1 = malloc(5);
            char *buf2 = malloc(5);
            if (!buf1 || !buf2) {
                perror("malloc");
                exit(EXIT_FAILURE);
            }
            buf1[0] = '\0';
            buf2[0] = '\0';
            
            strncat(buf1, source, 4);      // primele 4 biți
            strncat(buf2, source + 4, 4);  // următorii 4 biți
            
            int row2 = BinaryStringToInt(buf1);
            int col2 = BinaryStringToInt(buf2);
            
            free(buf1);
            free(buf2);
            
            uint8_t val = sbox_matrix[row2][col2];
            binary_state[k][j] = decToBinary(val);
        }
    }
}

void ShiftRows(char*** binary_state) {
    char *temp = NULL;
    
    printf("After ShiftRows:\n");
    
    // Rotate first row 1 column to left
    temp = binary_state[1][0];
    binary_state[1][0] = binary_state[1][1];
    binary_state[1][1] = binary_state[1][2];
    binary_state[1][2] = binary_state[1][3];
    binary_state[1][3] = temp;
    
    // Rotate second row 2 columns to left
    temp = binary_state[2][0];
    binary_state[2][0] = binary_state[2][2];
    binary_state[2][2] = temp;
    
    temp = binary_state[2][1];
    binary_state[2][1] = binary_state[2][3];
    binary_state[2][3] = temp;
    
    // Rotate third row 3 columns to left
    temp = binary_state[3][0];
    binary_state[3][0] = binary_state[3][3];
    binary_state[3][3] = binary_state[3][2];
    binary_state[3][2] = binary_state[3][1];
    binary_state[3][1] = temp;
}

int main(){
    int i=0, j=0, k=0;

    char *buf1;
    const char *message = "aaaabbbbccccdddd";
    const char *key = "cheie pe 16 biti"; 
    size_t size_message = strlen(message);
    size_t size_key_message = strlen(key);
    
    //state matrix is 4x4
    int rows = 4;
    int cols = 4;
    
    char** state = (char**)malloc(rows * sizeof(char*));
    for (i = 0; i < rows; i++) {
        state[i] = (char*)malloc(cols * sizeof(char));
    }

    char** w = (char**)malloc(rows * sizeof(char*));
    for (i = 0; i < rows; i++) {
        w[i] = (char*)malloc(cols * sizeof(char));
    }

    int counter = 0;

    // Create and initialize matrix (column-major order)
    for(j = 0; j < cols; j++) {
        for(i = 0; i < rows; i++) {
            if (counter < size_message){
                state[i][j] = message[counter];
                counter++;
            }
            else
                state[i][j] = 0x0F;  // padding
        }
    }

    //binary state
    char*** binary_state = (char***)malloc(rows * sizeof(char**));
    for (i = 0; i < rows; i++) {
        binary_state[i] = (char**)malloc(cols * sizeof(char*));
    }

    printf("binary_state: ");
    printf("\n");
    for(i = 0; i < rows; i++) {
        for(j = 0; j < cols; j++){
            int dec_value = (unsigned char)state[i][j];
            binary_state[i][j] = decToBinary(dec_value);
            printf("%s ", binary_state[i][j]);
        }
        printf("\n");
    }

    counter = 0;
    // Create and initialize encryption key matrix
    for(j = 0; j < cols; j++) {
        for(i = 0; i < rows; i++) {
            if (counter < size_key_message){
                w[i][j] = key[counter];
                counter++;
            }
        }
    }

    //binary state key
    char*** binary_state_key = (char***)malloc(rows * sizeof(char**));
    for (i = 0; i < rows; i++) {
        binary_state_key[i] = (char**)malloc(cols * sizeof(char*));
    }

    printf("binary_state_key: ");
    printf("\n");
    for(i = 0; i < rows; i++) {
        for(j = 0; j < cols; j++){
            int dec_value = (unsigned char)w[i][j];
            binary_state_key[i][j] = decToBinary(dec_value);
            printf("%s ", binary_state_key[i][j] );
        }
        printf("\n");
    }

    int expanded_cols = 8; // For key expansion
    char*** w2 = (char***)malloc(rows * sizeof(char**));
    for (i = 0; i < rows; i++) {
        w2[i] = (char**)malloc(expanded_cols * sizeof(char*));
        for (j = 0; j < expanded_cols; j++) {
            w2[i][j] = NULL;
        }
    }

    for(j = 0; j < 4; j++) {
        for(i = 0; i < rows; i++) {
            w2[i][j] = binary_state_key[i][j];
        }
    }

    printf("expansion key w2[0..3]: ");
    printf("\n");
    for(i = 0; i < 4; i++) {
        for(j = 0; j < 4; j++){
            printf("%s ", w2[i][j]);
        }
        printf("\n");
    }

    uint8_t sbox[256];
    initialize_aes_sbox(sbox);

    uint8_t sbox_matrix[16][16];
	j = 0;
	k = 0;
    for(i = 0; i < 256; i++){
		sbox_matrix[j][k] = sbox[i];
		k++;
		if ((i + 1) % 16 == 0){
			j++;
            k=0;
        }
    }

	//for(int i=0; i<16; i++){
	//	for(j=0; j<16; j++){
	//		printf("%02x ", sbox_matrix[i][j]);
	//	}
	//	printf("\n");
	//}

    AddRoundKey(binary_state, w2, sbox_matrix, expanded_cols, rows);

    printf("column 7 of w2 matrix is w_8");
    printf("\n");
    for(i = 0; i < 4; i++) {
        for(j = 4; j < expanded_cols; j++){
            if (w2[i][j])
                printf("%s ", w2[i][j]);
            else
                printf("NULL ");
        }
        printf("\n");
    }

    printf("the new binary_state after AddRoundKey: ");
    printf("\n");
    for(i = 0; i < 4; i++) {
        for(j = 0; j < 4; j++){
                printf("%s ", binary_state[i][j]);
        }
        printf("\n");
    }

    printf("Now begins the rounds from 1 to Nr-1 do: ");

    for(int round=1; round<Nr; round++){
        printf("\n=== Round %d ===\n", round);

        //S-box substitution
        SubBytes(binary_state, sbox_matrix);

        //SHIFTROWS()
        ShiftRows(binary_state);

        //MIXCOLUMNS
        char* mix_poly_binary[4][4] = {
            {"00000010", "00000011", "00000001", "00000001"},
            {"00000001", "00000010", "00000011", "00000001"},
            {"00000001", "00000001", "00000010", "00000011"},
            {"00000011", "00000001", "00000001", "00000010"}
        };

        char* new_binary_state[4][4];
        
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                char* rezultat_final[4] = {NULL, NULL, NULL, NULL};

                for (k = 0; k < 4; k++) {
                    rezultat_final[k] = multiplication_in_GF_2_power_8(mix_poly_binary[i][k], binary_state[k][j]); 
                }

                char* xor_temp = rezultat_final[0];

                for (k = 0; k < 4; k++) {
                    char* temp = xoring(xor_temp, rezultat_final[k], 8);
                    if (k > 1) 
                        free(xor_temp);
                    xor_temp = temp;
                }
                new_binary_state[i][j] = xor_temp;
                
                for (k = 0; k < 4; k++) {
                    if (k != 0) free(rezultat_final[k]);
                }
            }
        }
        printf("the new_binary_state after MixColumns: \n");
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                binary_state[i][j] = new_binary_state[i][j];
                printf("%s\t", binary_state[i][j]);
            }
            printf("\n");
        }

        //mai intai ne asiguram ca w3 are size [0..7]
        char*** w3 = (char***)malloc(rows * sizeof(char**));
        for (int it1 = 0; it1 < rows; it1++) {
            w3[it1] = (char**)malloc(expanded_cols * sizeof(char*));
            for (int it2 = 0; it2 < expanded_cols; it2++) {
                w3[it1][it2] = NULL;
            }
        }

        //acum copiem efectiv in w3[0..3] pe w2[4..7]
        for (int it1 = 0; it1 < rows; it1++) {
            for (int it2 = 4; it2 < expanded_cols; it2++) {
                w3[it1][it2-4] = w2[it1][it2];
            }
        }

        AddRoundKey(binary_state, w3, sbox_matrix, expanded_cols, rows);

        //actualizam w2 pentru urmatoarea iteratie
        for (int it1 = 0; it1 < rows; it1++) {
            for (int it2 = 0; it2 < expanded_cols; it2++) {
                w2[it1][it2] = w3[it1][it2];
            }
        }

        free(w3);
    }

    printf("the new binary_state after rounds from 1 to Nr-1: \n");
    printf("\n");
    for(i = 0; i < 4; i++) {
        for(j = 0; j < 4; j++){
            printf("%s ", binary_state[i][j]);
        }
        printf("\n");
    }

    printf("\n=== Round %d ===\n", 10);

    //S-box substitution
    SubBytes(binary_state, sbox_matrix);

    //SHIFTROWS()
    ShiftRows(binary_state);

    //mai intai ne asiguram ca w3 are size [0..7]
    char*** w3 = (char***)malloc(rows * sizeof(char**));
    for (int it1 = 0; it1 < rows; it1++) {
        w3[it1] = (char**)malloc(expanded_cols * sizeof(char*));
        for (int it2 = 0; it2 < expanded_cols; it2++) {
            w3[it1][it2] = NULL;
        }
    }

    //acum copiem efectiv in w3[0..3] pe w2[4..7]
    for (int it1 = 0; it1 < rows; it1++) {
        for (int it2 = 4; it2 < expanded_cols; it2++) {
            w3[it1][it2-4] = w2[it1][it2];
        }
    }

    AddRoundKey(binary_state, w3, sbox_matrix, expanded_cols, rows);

    //actualizam w2 pentru urmatoarea iteratie
    for (int it1 = 0; it1 < rows; it1++) {
        for (int it2 = 0; it2 < expanded_cols; it2++) {
            w2[it1][it2] = w3[it1][it2];
        }
    }

    free(w3);

    printf("binary_state dupa ce am facut de Nr+1 ori AddRoundKey: \n");
    printf("\n");
    for(i = 0; i < 4; i++) {
        for(j = 0; j < 4; j++){
            printf("%s ", binary_state[i][j]);
        }
        printf("\n");
    }

    // Free the memory
    for(i = 0; i < rows; i++) {
        free(state[i]);
    }

    for (i = 0; i < rows; i++) {
        for (j = 0; j < expanded_cols; j++) {
            free(w2[i][j]);
        }
        free(w2[i]);
    }
    free(w2);

    return 0;
}

