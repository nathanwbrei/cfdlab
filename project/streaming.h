#ifndef _STREAMING_H_
#define _STREAMING_H_

/** carries out the streaming step and writes the respective distribution functions from
 *  collideField to streamField.
 */
void doStreaming(float *collideField, float *streamField, int *flagField, float * massField, float * fractionField, int * length, int n_threads, float exchange);

#endif

