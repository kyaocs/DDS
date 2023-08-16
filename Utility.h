//
//  Utility.h
//  directed_densest_subgraph
//
//  Created by kai on 2023/8/7.
//

#ifndef Utility_h
#define Utility_h

#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <vector>
#include <queue>
//#define NDEBUG
#include <assert.h>
#include <algorithm>
#include <string.h>
#include <iomanip>
#include <math.h>
#include "Timer.h"

using namespace std;
typedef unsigned int ui;
const int INF = 1000000000;
#define pb push_back
#define mp make_pair

//#define _DEBUG_

ui n, m, n1;
ui * pstart;
ui * edges;
ui * degree;
bool * L;
ui * process_order;
ui max_core;
ui * core;
ui L_dmax, R_dmax;

long long bestL = 0, bestR = 0, bestM = 0;

string integer_to_string(long long number) {
    std::vector<ui> sequence;
    if(number == 0) sequence.push_back(0);
    while(number > 0) {
        sequence.push_back(number%1000);
        number /= 1000;
    }
    
    char buf[5];
    std::string res;
    for(unsigned int i = (ui)sequence.size();i > 0;i --) {
        if(i == sequence.size()) sprintf(buf, "%u", sequence[i-1]);
        else sprintf(buf, ",%03u", sequence[i-1]);
        res += std::string(buf);
    }
    return res;
}

#endif /* Utility_h */
