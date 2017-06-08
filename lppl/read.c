/** \brief I/O utilities

           Read price vector, safely open files */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "shared.h"


/** Open a file with error handling */
FILE *safe_fopen(const char *filename, char *error_msg) {
  FILE *fileptr;

  fileptr = fopen(filename, "r");
  perror_fail(fileptr == NULL, error_msg);

  return fileptr;
}


/** Open a command pipe with error handling */
FILE *safe_popen(char *command, char *type, char *error_msg) {
  FILE *fileptr;

  fileptr = popen(command, type);
  perror_fail(fileptr == NULL, error_msg);

  return fileptr;
}


/** Close a command pipe with error handling */
void safe_pclose(FILE *fileptr, char *error_msg) {
  if (pclose(fileptr) < 0) 
    perror(error_msg);   
}


/** Read prices from the pricefile. Prices are expected in the format
    index price
    @pre pricefile is not NULL
    @return Number of entries read */
unsigned int read_price(unsigned int n, /**< Price vector size */
                        double price[], /**< Price vector */
                        FILE *pricefile /**< File with the price vector */) {
  unsigned int i, j;

  assert(pricefile != NULL);

  for (i = 0; 
       i < n && fscanf(pricefile, "%u %lg", &j, price+i) == 2;
       i++);

  return i;
} 
