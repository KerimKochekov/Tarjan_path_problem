#include "template.h"
#include "RegEx.h"

using namespace std;

class Tarjan{
    private:
    int n, m, root;
    vector<vector<pair<int,int> > > out, in;
    vector<vector<int> > cycle, noncycle;
    vector<bool> removed;
    vector<int> order, header, ancestor;
    vector<int> h, e;
    vector<RegEx*> S;
    // vector<vector<RegEx*> > P;
    vector<RegEx*> P;
    
    void print_out_graph(){
        for (int i = 1; i <= n; i++){
            for (auto to: out[i])
                cout << i << " "<< to.ff << endl;
        }
        cout<<"-------------------(out)-----------------"<<endl;
    }

    void print_in_graph(){
        for (auto x: in[4])
            cout<<x.ff<<","<<x.ss<<endl;
        cout<<"---------"<<endl;
        for (int i = 1; i <= n; i++){
            for (auto to: in[i])
                cout << i << " "<< to.ff << " "<<to.ss<<endl;
        }
        cout<<"-------------------(in)-----------------"<<endl;
    }
    

    public:
    Tarjan(int _n, int _m){
        n = _n;
        m = _m;
        out.assign(n + 1, vector<PII>());
        in.assign(n + 1, vector<PII>());
        cycle.assign(n + 1, vector<int>());
        noncycle.assign(n + 1, vector<int>());
        removed.assign(n + 1, false);
        header.assign(n + 1, 0);
        ancestor.assign(n + 1, 0);
        S.assign(n + 1, NULL);
        order.clear();
        h.assign(m + 1, 0);
        e.assign(m + 1, 0);
        // P.assign(n + 1, vector<RegEx*> ());
        P.assign(n + 1, NULL);
    }

    void add_edge(int u, int v, int i){
        out[u].pb(mp(v, i));
        in[v].pb(mp(u, i));
        h[i] = u;
        e[i] = v;
    }

    void remove_out_edge(int v, pair<int, int> x){
        for (int i = 0; i< sz(out[v]); i++)
            if (out[v][i] == x){
                out[v].erase(out[v].begin() + i);
                break;
            }
    }

    void remove_in_edge(int v, pair<int, int> x){
        for (int i = 0; i< sz(in[v]); i++)
            if (in[v][i] == x){
                in[v].erase(in[v].begin() + i);
                break;
            }
    }

    void detach_node(int v){
        out[v].clear();
        in[v].clear();
    }

    // Algorithm starts here
    void update(int w, int v, RegEx *R){
        ancestor[v] = w;
        S[v] = R;
    }

    pair<RegEx*, int> eval(int v){
        if (ancestor[v] == 0)
            return mp(new RegEx(RegEx::ONE), v);
        pair<RegEx*, int> parent = eval(ancestor[v]);
        ancestor[v] = parent.ss;
        S[v] = new RegEx(RegEx::DOT, parent.ff, S[v]);
        return mp(S[v], ancestor[v]);
    }

    void reduce(){
        // initialize
        for (int v = 1; v <= n; v++){
            ancestor[v] = 0;
            S[v] = new RegEx(RegEx::ONE);
        }
        // loop
        for (int i = 0; i < n - 1; i++){
            int v = order[i];
            RegEx* P = new RegEx(RegEx::ZERO);
            RegEx* Q = new RegEx(RegEx::ZERO);
            for (auto e: noncycle[v])
                P = new RegEx(RegEx::PLUS, P, 
                    new RegEx(RegEx::DOT, eval(h[e]).ff, new RegEx(e)));
            for (auto e: cycle[v])
                Q = new RegEx(RegEx::PLUS, Q, 
                    new RegEx(RegEx::DOT, eval(h[e]).ff, new RegEx(e)));
            update(header[v], v, new RegEx(RegEx::DOT, P, new RegEx(RegEx::STAR, Q, nullptr)));
        }
        // finalize
        P[root] = new RegEx(RegEx::ZERO);
        for (auto e: cycle[root])
            P[root] = new RegEx(RegEx::PLUS, P[root], 
                    new RegEx(RegEx::DOT, eval(h[e]).ff, new RegEx(e)));
        P[root] = new RegEx(RegEx::STAR, P[root], nullptr);
        for (int i = 0; i < n - 1; i++)
            P[order[i]] = new RegEx(RegEx::DOT, P[root], eval(order[i]).ff);
    }

    void error_verdict(){
        puts("Graph is not reducible or flow graph");
        assert(0);
    }
    
    void prep(){
        for (int step = 1; step < n; step++){
            //T2
            int W = -1, V = - 1;
            for (int w = 1; w <= n; w++)
                if (!removed[w]){
                    if (in[w].empty())
                        error_verdict();
                    int v = in[w][0].ff;
                    bool t2 = true;
                    for (auto to: in[w])
                        t2 &= (to.ff == v);
                    if(t2){
                        W = w;
                        V = v;
                        break;
                    }
                }
            if (W == -1)
                error_verdict();
            order.pb(W);
            //T2 contraction
            for (auto to: in[W]){
                remove_out_edge(V, mp(W, to.ss));
                noncycle[W].pb(to.ss);
            }
            for (auto to: out[W]){
                out[V].pb(to);
                in[to.ff].pb(mp(V, to.ss));
                remove_in_edge(to.ff, mp(W, to.ss));
            }
            detach_node(W);
            header[W] = V;
            removed[W] = true;
            // print_in_graph();
            // print_out_graph();
            //T1
            bool go;
            do{
                go = false;
                for (auto to: out[V])
                    if (to.ff == V){
                        cycle[V].pb(to.ss);
                        remove_out_edge(V, mp(V, to.ss));
                        remove_in_edge(V, mp(V, to.ss));
                        go = true;
                        break;
                    }
            }while(go);
        }
        for (int i = 1; i <= n; i++)
            if (!removed[i])
                root = i;
        order.pb(root);
    }

    void print(){
        assert(order.size() == n);
        printf("order = [");
        for (int i = 0; i < n; i++)
            printf("%d,", order[i]);
        puts("]\n");
        for (int v = 1; v <= n; v++){
            printf("header(%d) = %d | ", v, header[v]);
            for (auto e: cycle[v])
                printf("e%d,", e);
            printf(" | ");
            for (auto e: noncycle[v])
                printf("e%d,", e);
            puts("");
        }
        puts("");
    }

    void print_regex(){
        for (int v = 1; v <= n; v++){
            if (v == root)
                cout << "root " << v << ": ";
            else
                cout << "node " << v << ": ";
            if (P[v] != NULL)
                cout << *P[v] << "\n";
            else    
                cout << "\n";
        }
    }
};
int main(){
    freopen("graph.txt", "r", stdin);
    int n, m;
    scanf("%d%d", &n, &m);
    Tarjan T(n, m);
    for (int i = 1; i <= m; i++){
        int u, v;
        scanf("%d%d", &u, &v);
        T.add_edge(u, v, i);
    }
    T.prep();
    T.print();
    T.reduce();
    T.print_regex();
    return 0;
}
