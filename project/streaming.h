#ifndef _STREAMING_H_
#define _STREAMING_H_

/** carries out the streaming step and writes the respective distribution functions from
 *  collideField to streamField.
 */
void doStreaming(   float *collideField,    /*Pointer to the start of the collide field, the data here will be read but not written*/
                    float *streamField,     /*Pointer to the start of the stream field, the data here will be written*/
                    int *flagField,         /*Pointer to the start of the flag field, the data here will be read but not written*/
                    float * massField,      /*Pointer to the start of the mass field, the data here will be read and written*/
                    float * fractionField,  /*Pointer to the start of the fraction field, the data here will be read but not written*/
                    int * length,           /*Interior length of the cavity */
                    int n_threads,          /*Number of threads to be used */
                    float exchange);        /*Exchange factor, this can make more splashing effects but can introduce instabilities*/

#endif
