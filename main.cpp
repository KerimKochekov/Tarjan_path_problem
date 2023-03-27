#include "template.h"
#include "RegEx.h"

using namespace std;

class PathSequence{
    private:
    RegEx* P;
    int u, v;
    public:
    PathSequence(RegEx* _P, int _u, int _v): 
        P(_P), u(_u), v(_v) {}
    void print(){
        cout << *P << " " << u << "-" << v << endl;
    }
    bool equal(int _u, int _v){
        return (u == _u and v == _v);
    }
    RegEx* getP(){
        return P;
    }
    int getU(){
        return u;
    }
    int getV(){
        return v;
    }
};

class Tarjan{
    private:
    int n, m, source, T;
    vector<vector<pair<int,int> > > out, in;
    vector<vector<int> > cycle, noncycle;
    vector<bool> removed;
    vector<int> order, decompose_order, header, ancestor;
    
    vector<vector<RegEx*> > p;
    vector<RegEx*> S;
    vector<PathSequence> sequence;
    vector<vector<PII> > adj;
    vector<vector<int> > anc, children, tree, nontree, sibling;
    vector<vector<int> > G, G_inv, SCC;
    vector<int> vis, idom, depth, h, t, id, sccOrder, sccId;
    
    void print_out_graph(){
        for (int i = 1; i <= n; i++){
            for (auto to: out[i])
                cout << i << " "<< to.ff << endl;
        }
        cout<<"-------------------(out)-----------------"<<endl;
    }

    void print_in_graph(){
        for (int i = 1; i <= n; i++){
            for (auto to: in[i])
                cout << i << " "<< to.ff << " "<<to.ss<<endl;
        }
        cout<<"-------------------(in)-----------------"<<endl;
    }
    
    void print_adj_graph(){
        for (int i = 1; i <= n; i++){
            for (auto to: adj[i])
                cout << i << " "<< to.ff << " "<<to.ss<<endl;
        }
        cout<<"-------------------(adj)-----------------"<<endl;
    }
    
    public:
    void initiate(int _n, int _m, int _source){
        n = _n;
        m = _m;
        source = _source;
        out.assign(n + 1, vector<PII>());
        in.assign(n + 1, vector<PII>());
        cycle.assign(n + 1, vector<int>());
        noncycle.assign(n + 1, vector<int>());
        removed.assign(n + 1, false);
        header.assign(n + 1, 0);
        ancestor.assign(n + 1, 0);
        S.assign(n + 1, NULL);
        sccId.assign(n + 1, 0);
        h.assign(m + 1, 0);
        t.assign(m + 1, 0);
        adj.assign(n + 1, vector<PII>());
        sibling.assign(n+1, vector<int>());
        anc.assign(n + 1, vector<int>());
        children.assign(n + 1, vector<int>());
        tree.assign(n + 1, vector<int>());
        nontree.assign(n + 1, vector<int>());
        G.assign(n + 1, vector<int>());
        G_inv.assign(n + 1, vector<int>());
        SCC.clear();
        idom.assign(n + 1, 0);
        depth.assign(n + 1, 0);
        id.assign(n + 1, 0);
    }

    vector<RegEx*> get_answer(int _n, int _m, int _source, vector<pair<int,int> > edges){
        initiate(_n, _m, _source);
        //make index one-based
        for (int i = 0; i < m; i++)
            add_edge(edges[i].ff, edges[i].ss, i + 1);
        confirm();
        compute();
        decompose_and_sequence();
        return solve(source);
    }

    void add_edge(int u, int v, int i){
        out[u].pb(mp(v, i));
        adj[u].pb(mp(v, i));
        in[v].pb(mp(u, i));
        h[i] = u; t[i] = v;
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
    vector<RegEx*> solve(int source){
        vector<RegEx*> Pt(n + 1, nullptr);
        Pt[source] = new RegEx(RegEx::ONE);
        for (int v = 1; v <= n; v++)
            if (v != source)
                Pt[v] = new RegEx(RegEx::ZERO);
        for (auto seq: sequence){
            int v = seq.getU();
            int w = seq.getV();
            if (v == w)
                Pt[v] = new RegEx(RegEx::DOT, Pt[v], seq.getP());
            else
                Pt[w] = new RegEx(RegEx::PLUS, Pt[w], new RegEx(RegEx::DOT, Pt[v], seq.getP()));
        }
        return Pt;
    }

    void show_answer(vector<RegEx*>p, string task_name){
        cout << "--------- " << task_name << " ---------" << endl;
        for (int v = 1; v <= n; v++){
            if (v == source)
                cout << "source " << v << ": ";
            else
                cout << "node " << v << ": ";
            if (p[v] != NULL)
                cout << *p[v] << "\n";
            else    
                cout << "\n";
        }
        puts("--------------------------");
    }

    vector<PathSequence> eleminate(vector<pair<PII,RegEx*> > edges){
        if (edges.empty())
            return {};
        map<int, int> id, who;
        int n_edges = edges.size(), vertices = 0;
        for (auto edge: edges){
            int u = edge.ff.ff;
            int v = edge.ff.ss;
            if (id.find(u) == id.end()){
                id[u] = vertices;
                who[vertices] = u;
                vertices += 1;
            }
            if (id.find(v) == id.end()){
                id[v] = vertices;
                who[vertices] = v;
                vertices += 1;
            }
        }
        p.assign(vertices, vector<RegEx*>(vertices, NULL));
        for (int v = 0; v < vertices; v++)
            for (int w = 0; w < vertices; w++)
                p[v][w] = new RegEx(RegEx::ZERO);
        for (auto edge: edges){
            int u = id[edge.ff.ff];
            int v = id[edge.ff.ss];
            p[u][v] = new RegEx(RegEx::PLUS, p[u][v], edge.ss);
        }
        for (int v = 0; v < vertices; v++){
            p[v][v] = new RegEx(RegEx::STAR, p[v][v], nullptr);
            for (int u = v + 1; u < vertices; u++){
                if (p[u][v] -> eId != RegEx::ZERO){
                    p[u][v] = new RegEx(RegEx::DOT, p[u][v], p[v][v]);
                    for (int w = v + 1; w < vertices; w++)
                        if (p[v][w] -> eId != RegEx::ZERO){
                            p[u][w] = new RegEx(RegEx::PLUS, p[u][w], 
                                new RegEx(RegEx::DOT, p[u][v], p[v][w]));
                        }
                }   
            }
        }
        vector<PathSequence> Y;
        for (int u = 0; u < vertices; u++)
            for (int v = u; v < vertices; v++)
                if (p[u][v] -> eId != RegEx::ZERO)
                    Y.pb(PathSequence(p[u][v], who[u], who[v]));
        for (int u = vertices - 1; u >= 0; u--)
            for (int v = 0; v < u; v++)
                if (p[u][v] -> eId != RegEx::ZERO)
                    Y.pb(PathSequence(p[u][v], who[u], who[v]));
        return Y;
    }
    // Runs eleminate function on entire graph
    void eleminate_test(bool show = false){
        vector<pair<PII,RegEx*> > edges;
        for (int i = 1; i <= m; i++)
            edges.pb(mp(mp(h[i],t[i]), get_edge(i)));
        sequence = eleminate(edges);
        if (show)
            show_answer(solve(source), "ELEMINATE");
    }

    void init(){
        for (int v = 1; v <= n; v++){
            ancestor[v] = 0;
            S[v] = new RegEx(RegEx::ONE);
        }
    }

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

    void reduce(bool show = false){
        init();
        // loop
        for (int i = 0; i < n - 1; i++){
            int v = order[i];
            RegEx* P = new RegEx(RegEx::ZERO);
            RegEx* Q = new RegEx(RegEx::ZERO);
            for (auto e: noncycle[v])
                P = new RegEx(RegEx::PLUS, P, 
                    new RegEx(RegEx::DOT, eval(h[e]).ff, get_edge(e)));
            for (auto e: cycle[v])
                Q = new RegEx(RegEx::PLUS, Q, 
                    new RegEx(RegEx::DOT, eval(h[e]).ff, get_edge(e)));
            update(header[v], v, new RegEx(RegEx::DOT, P, new RegEx(RegEx::STAR, Q, nullptr)));
        }
        // finalize
        vector<RegEx*> P(n+1, nullptr);
        P[source] = new RegEx(RegEx::ZERO);
        for (auto e: cycle[source])
            P[source] = new RegEx(RegEx::PLUS, P[source], 
                    new RegEx(RegEx::DOT, eval(h[e]).ff, get_edge(e)));
        P[source] = new RegEx(RegEx::STAR, P[source], nullptr);
        for (int i = 0; i < n - 1; i++)
            P[order[i]] = new RegEx(RegEx::DOT, P[source], eval(order[i]).ff);
        if (show)
            show_answer(P, "REDUCE");
    }
    
    pair<RegEx*, int> eval_and_compress(int v, int e){
        if (ancestor[v] == 0)
            return mp(new RegEx(RegEx::ONE), v);
        pair<RegEx*, int> parent = eval_and_compress(ancestor[v], e);
        ancestor[v] = parent.ss;
        S[v] = new RegEx(RegEx::DOT, parent.ff, S[v]);
        sequence.push_back(PathSequence(new RegEx(RegEx::DOT, S[v], get_edge(e)), v, t[e]));
        return mp(S[v], ancestor[v]);
    }

    RegEx* eval_and_sequence(int e){
        return new RegEx(RegEx::DOT, eval_and_compress(h[e], e).ff, get_edge(e));
    }

    void reduce_and_sequence(bool show = false){
        init();
        sequence.clear();
        for (int i = 0; i < n - 1; i++){
            int v = order[i];
            RegEx* P = new RegEx(RegEx::ZERO);
            RegEx* Q = new RegEx(RegEx::ZERO);
            for (auto e: noncycle[v])
                P = new RegEx(RegEx::PLUS, P, eval_and_sequence(e));
            for (auto e: cycle[v])
                Q = new RegEx(RegEx::PLUS, Q, eval_and_sequence(e));
            if (Q -> eId != RegEx::ONE and Q -> eId != RegEx::ZERO)
                sequence.pb(PathSequence(new RegEx(RegEx::STAR, Q, nullptr), v, v));
            update(header[v], v, new RegEx(RegEx::DOT, P, new RegEx(RegEx::STAR, Q, nullptr)));
        }
        RegEx* Q = new RegEx(RegEx::ZERO);
        for (auto e: cycle[source])
            Q = new RegEx(RegEx::PLUS, Q, eval_and_sequence(e));
        if (Q -> eId != RegEx::ONE and Q -> eId != RegEx::ZERO)
            sequence.pb(PathSequence(new RegEx(RegEx::STAR, Q, nullptr), source, source));
        
        for (int i = n - 2; i >= 0; i--){
            int v = order[i];
            sequence.pb(PathSequence(S[v], ancestor[v], v));
        }
        if (show)
            show_answer(solve(source), "REDUCE AND SEQUENCE");
    }
    
    void prep(){
        order.clear();
        for (int step = 1; step < n; step++){
            //T2
            int W = -1, V = - 1;
            for (int w = 1; w <= n; w++)
                if (!removed[w]){
                    assert (!in[w].empty());
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
            assert(~W);
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
                //T1 contraction
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
        order.pb(source);
    }

    //Dominator tree related functions starts here
    void reachable(int nd){
        if (vis[nd]) return;
        vis[nd] = 1;
        for (auto e: adj[nd])
            reachable(e.ff);
    }

    void confirm(){
        vis.assign(n + 1, 0);
        reachable(source);
        for (int i = 1; i <= n; i++)
            if (!vis[i]){
                puts("Your source node is invalid, choose different source to achieve flow graph");
                exit(1);
            }
    }

    void give_id(int nd){
        id[nd] = T--;
        for (auto to: children[nd])
            give_id(to);
    }
    
    bool is_ancestor(int u, int v){
        if (u == v)
            return true;
        for (auto x: anc[v])
            if (x == u)
                return true;
        return false;
    }

    //Dominator tree build in O(N*(N+M))
    void compute(){
        for (int v = 1; v <= n; v++)
            if (source != v){
                vis.assign(n + 1, 0);
                vis[v] = 1;
                reachable(source);
                for (int i = 1; i <= n; i++)
                    if (!vis[i]){
                        anc[i].pb(v);
                        depth[i] += 1;
                    }
            }
        for (int v = 1; v <= n; v++)
            if (v != source){
                int mx = -1;
                idom[v] = source;
                for (auto nd: anc[v])
                    if (depth[nd] > mx)
                        mx = depth[nd], idom[v] = nd;
            }
        assert(idom[source] == 0);
        for(int v = 1; v <= n; v++)
            children[idom[v]].pb(v);
        for (int e = 1; e <= m; e++){
            int u = t[e];
            if (h[e] == idom[u])
                tree[u].pb(e);
            else
                nontree[u].pb(e);
        }
        T = n;
        give_id(0);

        //get derived graph
        for (int e = 1; e <= m; e++){
            if (t[e] == source)
                sibling[source].pb(e);
            else if (h[e] != idom[t[e]]){
                for (auto u: children[idom[t[e]]])
                    if (is_ancestor(u, h[e]))
                        sibling[u].pb(e);
            }
        }

        //compute SCC for each sibling set
        vis.assign(n + 1, 0);
        vector<int> comp;
        for (int u = 0; u <= n; u++){
            sccOrder.clear();
            for (auto v: children[u])
                for (auto e: sibling[v]){
                    G[v].pb(t[e]);
                    G_inv[t[e]].pb(v);
                }
            for (auto v: children[u])
                if (!vis[v])
                    dfs1(v);
            for (auto v: children[u])
                vis[v] = 0;
            reverse(all(sccOrder));
            for (auto v: sccOrder)
                if (!vis[v]){
                    comp.clear();
                    dfs2(v, comp);
                    for (auto r: comp)
                        sccId[r] = SCC.size();
                    SCC.pb(comp);
                }
            for (auto v: children[u]){
                vis[v] = 0;
                G[v].clear();
                G_inv[v].clear();
            }
        }
    }

    //SCC dfs's
    void dfs1(int u){
        vis[u] = true;
        for (auto v: G[u])
            if (!vis[v])
                dfs1(v);
        sccOrder.pb(u);
    }
    
    void dfs2(int u, vector<int>&comp){
        vis[u] = true;
        comp.pb(u);
        for (auto v: G_inv[u])
            if (!vis[v])
                dfs2(v, comp);
    }
    
    RegEx* get_edge(int e){
        return new RegEx((((ll)(h[e])) << 32) | t[e]);
    }

    bool is_reducible(){
        for (auto scc: SCC)
            if (scc.size() > 1)
                return 0;
        return 1;
    }

    void decompose_and_sequence(bool show = false){
        init();
        sequence.clear();
        vector<pair<PII,RegEx*> >edges;
        vector<RegEx*> P(m + 1), R(n + 1);
        vector<int> decompose_order;
        for (int i = 1; i <= n; i++)
            decompose_order.pb(i);
        auto cmp = [&] (int x, int y){
            return id[x] < id[y];
        };
        sort(all(decompose_order), cmp);
        for (auto u: decompose_order){
            if (children[u].empty())
                continue;
            edges.clear();
            for (auto v: children[u])
                for (auto e: nontree[v])
                    P[e] = eval_and_sequence(e);
            for (auto v: children[u])
                for (auto e: sibling[v])
                    edges.pb(mp(mp(v, t[e]), P[e])); 
            vector<PathSequence> Yu = eleminate(edges);
            for (auto path: Yu)
                sequence.pb(path);
            for (auto v: children[u]){
                R[v] = new RegEx(RegEx::ZERO);
                for (auto e: tree[v])
                    R[v] = new RegEx(RegEx::PLUS, R[v], get_edge(e)); 
            }
            for (auto path: Yu){
                int w = path.getU();
                int x = path.getV();
                if (w == x)
                    R[w] = new RegEx(RegEx::DOT, R[w], path.getP());
                else
                    R[x] = new RegEx(RegEx::PLUS, R[x], new RegEx(RegEx::DOT, R[w], path.getP()));
            }
            for (auto v: children[u])
                update(u, v, R[v]);
        }
        RegEx* Q = new RegEx(RegEx::ZERO);
        for (auto e: nontree[source])
            Q = new RegEx(RegEx::PLUS, Q, eval_and_sequence(e));
        if (Q -> eId != RegEx::ONE and Q -> eId != RegEx::ZERO)
            sequence.pb(PathSequence(new RegEx(RegEx::STAR, Q, nullptr), source, source));
        for (int i = n - 2; i >= 0; i--){
            int v = decompose_order[i];
            sequence.pb(PathSequence(S[v], ancestor[v], v));
        }
        if (show)
            show_answer(solve(source), "DECOMPOSE AND SEQUENCE");
    }

    //Print part
    void show_auxiliary(){
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
};


int main(){
    freopen("graph.txt", "r", stdin);
    int n, m, r;
    scanf("%d%d%d", &n, &m, &r);
    assert(1 <= n);
    assert(1 <= r and r <= n);
    vector<PII> edges;
    for (int i = 1; i <= m; i++){
        int u, v;
        scanf("%d%d", &u, &v);
        edges.pb(mp(u, v));
    }
    Tarjan T;
    T.get_answer(n, m , r, edges);
    exit(0);
    //Check the given graph is flow graph from given source
    T.confirm();

    /* Naive O(n^3) algorithm test */
    T.eleminate_test();

    /* 
        Compute the Dominator tree along with SCC on siblings sets in O(NM) time
        It can be fastened if needed: https://dl.acm.org/doi/pdf/10.1145/357062.357071
    */
    T.compute();
    
    //If G is reducible with T1 and T2 contractions
    if (T.is_reducible()){
        // Pre-processing for reduce algorithm
        T.prep();
        // O(m log n) without path sequences
        T.reduce();
        // O(m logn) with path sequences
        T.reduce_and_sequence();
    }

    /*
        5. Decomposition using Dominators in O(m log n + t),
        where t is the time to find path sequences for the dominator strong components of G.
    */
    T.decompose_and_sequence();
}