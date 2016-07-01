#ifndef _STREAMING_H_
#define _STREAMING_H_

/** carries out the streaming step and writes the respective distribution functions from
 *  collideField to streamField.
 */
void doStreaming(float *collideField, float *streamField, int *flagField, double * massField, double * fractionField, int * length);

#endif

