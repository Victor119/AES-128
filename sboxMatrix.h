#pragma once

#include <stdint.h>
#include <stdio.h>

void initialize_aes_sbox(uint8_t sbox[256]);

// function to convert Hexadecimal to Binary Number
void HexToBin(char* hexdec);