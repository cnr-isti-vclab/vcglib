#ifndef NXS_REPORT_H
#define NXS_REPORT_H

#include <iostream>
using namespace std;

#include "stopwatch.h"

class Report {
public:
    void Start(unsigned int tot, unsigned int st) {
        total = tot;
        step = st;
        watch.Restart();
    }
    void Output(unsigned int count) {
        if((count%step) != 0) return;
        float r = 100.0 * count/total;
        float t = watch.Elapsed();
        cout << r << "% in " << t << "s -> " << count/t << "n/s\n";
    }
    double Elapsed() {
        return watch.Elapsed();
    }
private:
    StopWatch watch;
    unsigned int total; //when finished
    unsigned int step; //check once every step
};

#endif
