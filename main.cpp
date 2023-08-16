//
//  main.cpp
//  directed_densest_subgraph
//
//  Created by kai on 2023/8/7.
//

#include "Utility.h"
#include "LinearHeap.h"

void load_graph(string input_graph)
{
    string buffer;
    ifstream input_file(input_graph, ios::in);
    if (!input_file.is_open()){cout << "cannot open file : "<<input_graph<<endl;exit(1);}
    else{
        input_file >> n >> m;
        //transfer this directed graph to a bipartite graph
//        L = new bool[2*n];
        ui u, v;
        vector<vector<ui>> G;
        G.resize(2*n);
        while (input_file >> u >> v){
            if(u == v || u >= n) exit(1);
            G[u].push_back(n+v);
            G[n+v].push_back(u);
        }
        input_file.close();
        if(G.size() != 2*n) exit(1);
        pstart = new ui[2*n+1];
        edges = new ui[2*m];
        degree = new ui[2*n];
        pstart[0] = 0;
        for(ui i = 0; i < 2*n; i++) {
            ui s_idx = pstart[i];
            for(auto &e : G[i]) edges[s_idx++] = e;
            pstart[i+1] = s_idx;
            degree[i] = (ui)G[i].size();
        }
        L_dmax = 0;
        R_dmax = 0;
        for(ui i = 0; i < n; i++) if(degree[i] > L_dmax) L_dmax = degree[i];
        for(ui i = n; i < 2*n; i++) if(degree[i] > R_dmax) R_dmax = degree[i];
        cout<<"G:"<<input_graph<<", n="<<n<<", m="<<m<<", L_dmax="<<L_dmax<<", R_dmax="<<R_dmax<<endl;
    }
    n1 = n;
    n = 2*n;
//    cout<<"bipartite graph inforamtion: "<<endl;
//    for(ui i = 0; i < 2*n; i++) {
//        cout<<"vertex "<<i<<" ("<<degree[i]<<"): ";
//        for(ui j = pstart[i]; j < pstart[i+1]; j++) cout<<edges[j]<<","; cout<<endl;
//    }
}

double HApprox()
{
    double bestGeoDen = 0.0;
    long long bestGeoDen_L = 0;
    long long bestGeoDen_R = 0;
    long long bestGeoDen_M = 0;
    
//    cout<<"** H **"<<endl;
    long long cur_n = n;
    long long cur_m = m;
    long long curL = n1;
    long long curR = n1;
//    int bestL = 0;
//    int bestR = 0;
//    int bestM = 0;
    double cur_den = (double)cur_m/cur_n;
    bestL = curL;
    bestR = curR;
    bestM = cur_m;
    long long opt_pos = 0;
    process_order = new ui[n];
    for(ui i = 0; i < n; i++) process_order[i] = i;
    max_core = 0;
    core = new ui[n];
    memset(core, 0, sizeof(ui)*n);
    ListLinearHeap *linear_heap = new ListLinearHeap(n, n-1);
    linear_heap->init(n, n-1, process_order, degree);
    for(ui i = 0; i < n; i ++) {
        ui u, key;
        linear_heap->pop_min(u, key);
        if(key > max_core) max_core = key;
        process_order[i] = u;
        core[u] = max_core;
        -- cur_n;
        assert(u >= 0&& u < n);
        if(u < n1) -- curL;
        else -- curR;
        for(ui j = pstart[u];j < pstart[u+1];j ++) if(core[edges[j]] == 0) {
            linear_heap->decrement(edges[j], 1);
            -- cur_m;
        }
        if((double)cur_m/cur_n >= cur_den) {
            cur_den = (double)cur_m/cur_n;
            opt_pos = i+1;
            bestL = curL;
            bestR = curR;
            bestM = cur_m;
            assert(curL + curR == cur_n);
        }
        if((double)cur_m/sqrt(curL*curR) > bestGeoDen) {
            bestGeoDen = (double)cur_m/(sqrt(curL*curR));
            bestGeoDen_L = curL;
            bestGeoDen_R = curR;
            bestGeoDen_M = cur_m;
        }
    }
    delete linear_heap;
    cur_den *= 2; //need to time 2 in the case of directed graphs
    cout<<"Harmonic Density = "<<cur_den<<", ("<<bestL<<", "<<bestR<<", "<<bestM<<")"<<endl;
//    cout<<"** bestGeoDen = "<<bestGeoDen<<", ("<<bestGeoDen_L<<", "<<bestGeoDen_R<<", "<<bestGeoDen_M<<")"<<endl;
    
#ifdef _DEBUG_
    //the results is process_order[opt_pos -> n];
    cout<<"the subgraph induced by :"; for(ui i = opt_pos; i < n; i++) cout<<process_order[i]<<","; cout<<endl;
    cout<<"process_order:"; for(ui i = 0; i < n; i++) cout<<process_order[i]<<","; cout<<endl;
    cout<<"core_number  :"; for(ui i = 0; i < n; i++) cout<<core[i]<<","; cout<<endl;
#endif
    return cur_den;
}

int compute_gamma()
{
    process_order = new ui[n];
    for(ui i = 0; i < n; i++) process_order[i] = i;
    max_core = 0;
    core = new ui[n];
    memset(core, 0, sizeof(ui)*n);
    ui * tmp_degree = new ui[n];
    for(ui i = 0; i < n; i++) tmp_degree[i] = degree[i];
    ListLinearHeap *linear_heap = new ListLinearHeap(n, n-1);
    linear_heap->init(n, n-1, process_order, tmp_degree);
    for(ui i = 0; i < n; i ++) {
        ui u, key;
        linear_heap->pop_min(u, key);
        if(key > max_core) max_core = key;
        process_order[i] = u;
        core[u] = max_core;
        assert(u >= 0&& u < n);
        for(ui j = pstart[u];j < pstart[u+1];j ++) if(core[edges[j]] == 0) {
            linear_heap->decrement(edges[j], 1);
        }
    }
    delete [] process_order;
    delete [] core;
    delete linear_heap;
    return max_core;
}

double GApprox()
{
//    cout<<"** G **"<<endl;
    int gamma = compute_gamma();
    cout<<"gamma = "<<gamma<<endl;
    
    int opt_tim = 1;
    int xstar = 0;
    int ystar = 0;
    int optL = 0;
    int optR = 0;
    int optM = 0;
//    int bestL = 0;
//    int bestR = 0;
//    int bestM = 0;
    
    ui * peel_order2;
    peel_order2 = new ui[n1];
    bool * del = new bool[n];
    ui * tmp_degree = new ui[n];
    ListLinearHeap *LH2 = new ListLinearHeap(n, R_dmax);
    //x
    for(int x = 1; x <= gamma; x++) {
//        cout<<"x="<<x<<endl;
        //prepare
        int cur_m = m;
        int ymax = 0;
        int Lsize = n1;
        int Rsize = n1;
        memset(del, 0, sizeof(bool)*n);
        for(ui i = 0; i < n; i++) tmp_degree[i] = degree[i];
        for(ui i = 0; i < n1; i++) peel_order2[i] = i + n1;
        LH2->init(n1, R_dmax, peel_order2, tmp_degree);
        
        //clear unqualified vertices in L
        for(ui i = 0; i < n1; i++) if(tmp_degree[i] < x) {
            del[i] = true;
            -- Lsize;
            for(ui j = pstart[i]; j < pstart[i+1]; j++) {
                LH2->decrement(edges[j], 1);
                -- cur_m;
            }
        }
        
        //find the largest y by increasing y from 1 to R_dmax
        for(int y = (opt_tim/x) + 1; y <= R_dmax + 1; y++) {
//            cout<<"\ty="<<y<<endl;
            ui tmp_v, tmp_deg;
            
            LH2->get_max(tmp_v, tmp_deg);
            if(tmp_deg < y) continue;
            
            LH2->get_min(tmp_v, tmp_deg);
            if(tmp_deg < y) {
                do {
                    LH2->pop_min(tmp_v, tmp_deg);
                    assert(tmp_deg < y);
                    del[tmp_v] = true;
                    -- Rsize;
                    for(ui i = pstart[tmp_v]; i < pstart[tmp_v+1]; i++) if(del[edges[i]] == false) {
                        assert(edges[i] < n1 && tmp_degree[edges[i]] > 0);
                        -- tmp_degree[edges[i]];
                        -- cur_m;
                        if(tmp_degree[edges[i]] < x) {
                            assert(tmp_degree[edges[i]] == x - 1);
                            del[edges[i]] = true;
                            -- Lsize;
                            for(ui j = pstart[edges[i]]; j < pstart[edges[i]+1]; j++) if(del[edges[j]] == false)
                            {
                                LH2->decrement(edges[j], 1);
                                -- cur_m;
                            }
                        }
                    }
                    LH2->get_min(tmp_v, tmp_deg);
                } while (tmp_deg < y && cur_m > 0);
            }
            if(cur_m > 0) {
                ymax = y; //exist such an (x,y)-core
                optM = cur_m;
                optL = Lsize;
                optR = Rsize;
            }
        }
//        assert(cur_m == 0);
        if(ymax*x > opt_tim) {
            opt_tim = ymax*x;
            xstar = x;
            ystar = ymax;
            bestM = optM;
            bestL = optL;
            bestR = optR;
        }
//        cout<<"opt_sum="<<opt_sum<<endl;
    }
    //y
    ui * peel_order1;
    peel_order1 = new ui[n1];
    ListLinearHeap *LH1 = new ListLinearHeap(n1, L_dmax);
    for(int y = 1; y <= gamma; y++) {
        //prepare
        int cur_m = m;
        int xmax = 0;
        int Lsize = n1;
        int Rsize = n1;
        memset(del, 0, sizeof(bool)*n);
        for(ui i = 0; i < n; i++) tmp_degree[i] = degree[i];
        for(ui i = 0; i < n1; i++) peel_order1[i] = i;
        LH1->init(n1, L_dmax, peel_order1, tmp_degree);
        
        //clear unqualified vertices in R
        for(ui i = n1; i < n; i++) if(tmp_degree[i] < y) {
            del[i] = true;
            -- Rsize;
            for(ui j = pstart[i]; j < pstart[i+1]; j++) {
                LH1->decrement(edges[j], 1);
                -- cur_m;
            }
        }
        
        //find the largest x by increasing x from 1 to L_dmax
        for(int x = (opt_tim/y) + 1; x <= L_dmax + 1; x++) {

            ui tmp_v, tmp_deg;
            
            LH1->get_max(tmp_v, tmp_deg);
            if(tmp_deg < x) continue;
            
            LH1->get_min(tmp_v, tmp_deg);
            if(tmp_deg < x) {
                do {
                    LH1->pop_min(tmp_v, tmp_deg);
                    assert(tmp_deg < x);
                    del[tmp_v] = true;
                    -- Lsize;
                    for(ui i = pstart[tmp_v]; i < pstart[tmp_v+1]; i++) if(del[edges[i]] == false) {
                        assert(edges[i] >= n1 && tmp_degree[edges[i]] > 0);
                        -- tmp_degree[edges[i]];
                        -- cur_m;
                        if(tmp_degree[edges[i]] < y) {
                            assert(tmp_degree[edges[i]] == y - 1);
                            del[edges[i]] = true;
                            -- Rsize;
                            for(ui j = pstart[edges[i]]; j < pstart[edges[i]+1]; j++) if(del[edges[j]] == false)
                            {
                                LH1->decrement(edges[j], 1);
                                -- cur_m;
                            }
                        }
                    }
                    LH1->get_min(tmp_v, tmp_deg);
                } while (tmp_deg < x && cur_m > 0);
            }
            if(cur_m > 0) {
                xmax = x; //exist such an (x,y)-core
                optM = cur_m;
                optL = Lsize;
                optR = Rsize;
            }
        }
        if(xmax*y > opt_tim) {
            opt_tim = xmax*y;
            xstar = xmax;
            ystar = y;
            bestM = optM;
            bestL = optL;
            bestR = optR;
        }
    }
    delete [] peel_order1;
    delete [] peel_order2;
    delete LH1;
    delete LH2;
    delete [] del;
    delete [] tmp_degree;
    double den = (double)bestM/sqrt(bestL*bestR);
    cout<<"Max-sum cn-pair: ("<<xstar<<","<<ystar<<")-core."<<endl;
    cout<<"Geometric Density = "<<den<<", ("<<bestL<<", "<<bestR<<", "<<bestM<<")"<<endl;
    return den;
}

double AApprox()
{
//    cout<<"** A **"<<endl;
    int gamma = compute_gamma();
    cout<<"gamma = "<<gamma<<endl;
    
    int opt_sum = 1;
    int xstar = 0;
    int ystar = 0;
    int optL = 0;
    int optR = 0;
    int optM = 0;
//    int bestL = 0;
//    int bestR = 0;
//    int bestM = 0;
    
    ui * peel_order2;
    peel_order2 = new ui[n1];
    bool * del = new bool[n];
    ui * tmp_degree = new ui[n];
    ListLinearHeap *LH2 = new ListLinearHeap(n, R_dmax);
    //x
    for(int x = 1; x <= gamma; x++) {
//        cout<<"x="<<x<<endl;
        //prepare
        int cur_m = m;
        int ymax = 0;
        int Lsize = n1;
        int Rsize = n1;
        memset(del, 0, sizeof(bool)*n);
        for(ui i = 0; i < n; i++) tmp_degree[i] = degree[i];
        for(ui i = 0; i < n1; i++) peel_order2[i] = i + n1;
        LH2->init(n1, R_dmax, peel_order2, tmp_degree);
        
        //clear unqualified vertices in L
        for(ui i = 0; i < n1; i++) if(tmp_degree[i] < x) {
            del[i] = true;
            -- Lsize;
            for(ui j = pstart[i]; j < pstart[i+1]; j++) {
                LH2->decrement(edges[j], 1);
                -- cur_m;
            }
        }
        
        //find the largest y by increasing y from 1 to R_dmax
        for(int y = opt_sum - x + 1; y <= R_dmax + 1; y++) {
//            cout<<"\ty="<<y<<endl;
            ui tmp_v, tmp_deg;
            
            LH2->get_max(tmp_v, tmp_deg);
            if(tmp_deg < y) continue;
            
            LH2->get_min(tmp_v, tmp_deg);
            if(tmp_deg < y) {
                do {
                    LH2->pop_min(tmp_v, tmp_deg);
                    assert(tmp_deg < y);
                    del[tmp_v] = true;
                    -- Rsize;
                    for(ui i = pstart[tmp_v]; i < pstart[tmp_v+1]; i++) if(del[edges[i]] == false) {
                        assert(edges[i] < n1 && tmp_degree[edges[i]] > 0);
                        -- tmp_degree[edges[i]];
                        -- cur_m;
                        if(tmp_degree[edges[i]] < x) {
                            assert(tmp_degree[edges[i]] == x - 1);
                            del[edges[i]] = true;
                            -- Lsize;
                            for(ui j = pstart[edges[i]]; j < pstart[edges[i]+1]; j++) if(del[edges[j]] == false)
                            {
                                LH2->decrement(edges[j], 1);
                                -- cur_m;
                            }
                        }
                    }
                    LH2->get_min(tmp_v, tmp_deg);
                } while (tmp_deg < y && cur_m > 0);
            }
            if(cur_m > 0) {
                ymax = y; //exist such an (x,y)-core
                optM = cur_m;
                optL = Lsize;
                optR = Rsize;
            }
        }
//        assert(cur_m == 0);
        if(ymax + x > opt_sum) {
            opt_sum = ymax + x;
            xstar = x;
            ystar = ymax;
            bestM = optM;
            bestL = optL;
            bestR = optR;
        }
//        cout<<"opt_sum="<<opt_sum<<endl;
    }
    //y
    ui * peel_order1;
    peel_order1 = new ui[n1];
    ListLinearHeap *LH1 = new ListLinearHeap(n1, L_dmax);
    for(int y = 1; y <= gamma; y++) {
        //prepare
        int cur_m = m;
        int xmax = 0;
        int Lsize = n1;
        int Rsize = n1;
        memset(del, 0, sizeof(bool)*n);
        for(ui i = 0; i < n; i++) tmp_degree[i] = degree[i];
        for(ui i = 0; i < n1; i++) peel_order1[i] = i;
        LH1->init(n1, L_dmax, peel_order1, tmp_degree);
        
        //clear unqualified vertices in R
        for(ui i = n1; i < n; i++) if(tmp_degree[i] < y) {
            del[i] = true;
            -- Rsize;
            for(ui j = pstart[i]; j < pstart[i+1]; j++) {
                LH1->decrement(edges[j], 1);
                -- cur_m;
            }
        }
        
        //find the largest x by increasing x from 1 to L_dmax
        for(int x = opt_sum - y + 1; x <= L_dmax + 1; x++) {

            ui tmp_v, tmp_deg;
            
            LH1->get_max(tmp_v, tmp_deg);
            if(tmp_deg < x) continue;
            
            LH1->get_min(tmp_v, tmp_deg);
            if(tmp_deg < x) {
                do {
                    LH1->pop_min(tmp_v, tmp_deg);
                    assert(tmp_deg < x);
                    del[tmp_v] = true;
                    -- Lsize;
                    for(ui i = pstart[tmp_v]; i < pstart[tmp_v+1]; i++) if(del[edges[i]] == false) {
                        assert(edges[i] >= n1 && tmp_degree[edges[i]] > 0);
                        -- tmp_degree[edges[i]];
                        -- cur_m;
                        if(tmp_degree[edges[i]] < y) {
                            assert(tmp_degree[edges[i]] == y - 1);
                            del[edges[i]] = true;
                            -- Rsize;
                            for(ui j = pstart[edges[i]]; j < pstart[edges[i]+1]; j++) if(del[edges[j]] == false)
                            {
                                LH1->decrement(edges[j], 1);
                                -- cur_m;
                            }
                        }
                    }
                    LH1->get_min(tmp_v, tmp_deg);
                } while (tmp_deg < x && cur_m > 0);
            }
            if(cur_m > 0) {
                xmax = x; //exist such an (x,y)-core
                optM = cur_m;
                optL = Lsize;
                optR = Rsize;
            }
        }
        if(xmax + y > opt_sum) {
            opt_sum = xmax + y;
            xstar = xmax;
            ystar = y;
            bestM = optM;
            bestL = optL;
            bestR = optR;
        }
    }
    delete [] peel_order1;
    delete [] peel_order2;
    delete LH1;
    delete LH2;
    delete [] del;
    delete [] tmp_degree;
    double den = (double)((bestL+bestR)*(bestM))/(2*bestL*bestR);
    cout<<"Max-sum cn-pair: ("<<xstar<<","<<ystar<<")-core."<<endl;
    cout<<"Arithmetic Density = "<<den<<", ("<<bestL<<", "<<bestR<<", "<<bestM<<")"<<endl;
    return den;
}

double MApprox()
{
    double bestGeoDen = 0.0;
    long long bestGeoDen_L = 0;
    long long bestGeoDen_R = 0;
    long long bestGeoDen_M = 0;
    
//    cout<<"** M **"<<endl;
    ui * peel_order1;
    peel_order1 = new ui[n1];
    ListLinearHeap *LH1 = new ListLinearHeap(n1, L_dmax);
    for(ui i = 0; i < n1; i++) peel_order1[i] = i;
    LH1->init(n1, L_dmax, peel_order1, degree);
    
    ui * peel_order2;
    peel_order2 = new ui[n1];
    ListLinearHeap *LH2 = new ListLinearHeap(n, R_dmax);
    for(ui i = 0; i < n1; i++) peel_order2[i] = i + n1;
    LH2->init(n1, R_dmax, peel_order2, degree);
    
    bool * del = new bool[n];
    memset(del, 0, sizeof(bool)*n);
    
    long long cur_n1 = n1;
    long long cur_m = m;
    double cur_den = (double)cur_m/cur_n1;
    long long bestN = 0;
//    int bestM = 0;
    
    for(ui i = 0; i < n1; i++) {
        ui u, udeg;
        LH1->pop_min(u, udeg);
        del[u] = true;
        for(ui j = pstart[u]; j < pstart[u+1]; j++) if(del[edges[j]] == false) {
            assert(edges[j] >= n1);
            LH2->decrement(edges[j], 1);
            -- cur_m;
        }
        ui v, vdeg;
        LH2->pop_min(v, vdeg);
        del[v] = true;
        for(ui j = pstart[v]; j < pstart[v+1]; j++) if(del[edges[j]] == false) {
            assert(edges[j] >= 0 && edges[j] < n1);
            LH1->decrement(edges[j], 1);
            -- cur_m;
        }
        -- cur_n1;
        if((double)cur_m/cur_n1 > cur_den) {
            cur_den = (double)cur_m/cur_n1;
            bestN = cur_n1;
            bestM = cur_m;
        }
        if((double)cur_m/(sqrt(cur_n1*cur_n1)) > bestGeoDen) {
            bestGeoDen = (double)cur_m/(sqrt(cur_n1*cur_n1));
            bestGeoDen_L = cur_n1;
            bestGeoDen_R = cur_n1;
            bestGeoDen_M = cur_m;
        }
    }
    assert(cur_n1 == 0);
    assert(cur_m == 0 || cur_m == 1);
    bestL = bestN;
    bestR = bestN;
    cout<<"Minimum Density = "<<cur_den<<", ("<<bestL<<", "<<bestR<<", "<<bestM<<")"<<endl;
//    cout<<"** bestGeoDen = "<<bestGeoDen<<", ("<<bestGeoDen_L<<", "<<bestGeoDen_R<<", "<<bestGeoDen_M<<")"<<endl;
    
    delete [] peel_order1;
    delete [] peel_order2;
    delete LH1;
    delete LH2;
    delete [] del;
    
    return cur_den;
}

int main(int argc, const char * argv[]) {
    
    load_graph(argv[1]);
    
    string algo = argv[2];
    
    Timer t;
    
    if(algo.compare("H") == 0) HApprox();
    
    else if (algo.compare("A") == 0) AApprox();
    
    else if (algo.compare("M") == 0) MApprox();
    
    else if (algo.compare("G") == 0) GApprox();
    
    else cout<<"no matched algorithm!"<<endl;
        
    cout<<"Time cost: "<<integer_to_string(t.elapsed())<<endl;
    
//    cout<<"H = "<<(double)(2*bestM)/(bestL+bestR)<<endl;
//    cout<<"A = "<<(double)((bestL+bestR)*(bestM))/(2*bestL*bestR)<<endl;
//    cout<<"M = "<<(double)(bestM)/(max(bestL, bestR))<<endl;
//    cout<<"G = "<<(double)bestM/sqrt(bestL*bestR)<<endl;
    
    return 0;
}
