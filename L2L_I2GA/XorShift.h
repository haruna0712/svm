#pragma once
#pragma once
//XorShift Module
#include <stdint.h>
#include <string>
#include <vector>
#include <ctime>

namespace XorShifts {
    class XorShift {
    private:
        struct dictionary { uint32_t x; uint32_t y; uint32_t z; uint32_t w; };
        uint32_t x;
        uint32_t y;
        uint32_t z;
        uint32_t w;

    public:
        static const struct dictionary defaults;
        uint32_t randCount = 0;
        struct dictionary seeds;

        XorShift(
            uint32_t w = time(nullptr),
            uint32_t x = NULL,
            uint32_t y = NULL,
            uint32_t z = NULL
        ) {
            x = x != NULL ? x : w << 13;
            y = y != NULL ? y : (w >> 9) ^ (x << 6);
            z = z != NULL ? z : y >> 7;
            seeds = { x,y,z,w };
            this->w = w; this->x = x; this->y = y; this->z = z;
        }

        uint32_t rand() {
            randCount++;
            uint32_t t = x ^ (x << 11);
            x = y;
            y = z;
            z = w;
            return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
        }

        double randDouble(double min = 0, double max = 1) {
            return (double)(rand() % 0xFFFF) / 0xFFFF * (max - min) + min;
        }


    };


}
//XorShift End