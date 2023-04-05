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
    int n, m, T, source;

    vector<int> out;
    vector<int> order, decompose_order, ancestor;
    vector<RegEx*> Pt;
    
    vector<vector<RegEx*> > p;
    vector<RegEx*> S;
    vector<PathSequence> sequence;
    vector<vector<PII> > adj;
    vector<vector<int> > anc, children, tree, nontree, sibling;
    vector<vector<int> > G, G_inv, SCC;
    vector<int> vis, idom, depth, h, t, reach;
    vector<int> topid, id, sccOrder, sccId, reach_edges;
    vector<vector<int> > D;

    vector<vector<pair<vector<pair<pair<int,int>,int>>,
                vector<pair<pair<int,int>,int>> > >> SC;
    
    void print_adj_graph(){
        for (int i = 1; i <= n; i++){
            for (auto to: adj[i])
                cout << i << " "<< to.ff << " "<<to.ss<<endl;
        }
        cout<<"-------------------(adj)-----------------"<<endl;
    }
    
    public:
    Tarjan(int _n, int _m, vector<pair<int,int> > edges){
        n = _n;
        m = _m;

        vis.assign(n + 1, 0);
        out.assign(n + 1, 0);
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
        topid.assign(n + 1, 0);
        reach.clear();
        reach_edges.clear();
        Pt.assign(n + 1, nullptr);
        D.assign(n + 1, vector<int>());
        SC.assign(n + 1, vector<pair<vector<pair<pair<int,int>,int>>,
                              vector<pair<pair<int,int>,int>> >>());
        
        for (int i = 0; i < m; i++)
            add_edge(edges[i].ff, edges[i].ss, i + 1);
    }

    RegEx* query(int u, int v){
        source = u;
        reachable(source);
        bool not_reachable = !vis[v];
        for (auto w: reach)
            vis[w] = 0;
        if (not_reachable)
            return new RegEx(RegEx::ZERO);
        compute();
        decompose_and_sequence();
        solve();
        // If you want to print answer for all reachable nodes
        // show_answer();
        return Pt[v];
    }

    void add_edge(int u, int v, int i){
        adj[u].pb(mp(v, i));
        h[i] = u; t[i] = v;
    }

    // Algorithm starts here
    void solve(){
        Pt[source] = new RegEx(RegEx::ONE);
        for (auto v: reach)
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
    }

    void show_answer(){
        for (auto v: reach){
            if (v == source)
                cout << "source " << v << ": ";
            else
                cout << "node " << v << ": ";
            if (Pt[v] != NULL)
                cout << *Pt[v] << "\n";
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
   
    void init(){
        for (auto v: reach){
            ancestor[v] = 0;
            S[v] = new RegEx(RegEx::ONE);
        }
    }

    void update(int w, int v, RegEx *R){
        ancestor[v] = w;
        S[v] = R;
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
    
    void reachable(int nd){
        vis[nd] = 1;
        reach.pb(nd);
        for (auto e: adj[nd]){
            reach_edges.pb(e.ss);
            if(!vis[e.ff])
                reachable(e.ff);
        }
    }

    void give_id(int nd){
        id[nd] = T--;
        for (auto to: children[nd])
            give_id(to);
    }
    
    void dfs0(int nd){
        vis[nd] = 1;
        for (auto e: adj[nd])
            if(!vis[e.ff])
                dfs0(e.ff);
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
    
    bool is_ancestor(int u, int v){
        if (u == v)
            return true;
        for (auto x: anc[v])
            if (x == u)
                return true;
        return false;
    }

    void compute(){
        for (auto v: reach)
            if (source != v){
                vis[v] = 1;
                dfs0(source);
                for (auto i: reach){
                    if (!vis[i]){
                        anc[i].pb(v);
                        depth[i] += 1;
                    }
                    else    
                        vis[i] = 0;
                }
            }
        for (auto v: reach)
            if (v != source){
                int mx = -1;
                idom[v] = source;
                for (auto nd: anc[v])
                    if (depth[nd] > mx)
                        mx = depth[nd], idom[v] = nd;
            }
        idom[source] = 0;
        for (auto v: reach)
            children[idom[v]].pb(v);
        for (auto e: reach_edges){
            int u = t[e];
            if (h[e] == idom[u])
                tree[u].pb(e);
            else
                nontree[u].pb(e);
        }
        T = n;
        give_id(0);

        //get derived graph
        for (auto e: reach_edges){
            if (t[e] == source)
                sibling[source].pb(e);
            else if (h[e] != idom[t[e]]){
                for (auto u: children[idom[t[e]]])
                    if (is_ancestor(u, h[e]))
                        sibling[u].pb(e);
            }
        }

        //compute SCC for each sibling set
        vector<int> comp;
        vector<vector<int> > SCCL, SCCR;
        queue<int> q;
        vector<pair<pair<int,int>, int> > X, Y;
        reach.pb(0);
        for (auto u: reach){
            sccOrder.clear();
            SCCL.clear();
            SCCR.clear();
            for (auto v: children[u])
                for (auto e: sibling[v]){
                    G[v].pb(e);
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
                        sccId[r] = SCCL.size();
                    SCCL.pb(comp);
                    SCC.pb(comp);
                }
            // Topological sorting on SCC transformed graph
            int nl = SCCL.size();
            for (auto v: children[u])
                for (auto e: sibling[v])
                    if (sccId[v] != sccId[t[e]]){
                        D[sccId[t[e]]].pb(sccId[v]);
                        out[sccId[v]] += 1;
                    }
            for (int i = 0; i < nl; i++)
                if (!out[i])
                    q.push(i);
            while (!q.empty()){
                int nd = q.front();
                q.pop();
                SCCR.pb(SCCL[nd]);
                for (auto to: D[nd])
                    if (!(--out[to]))
                        q.push(to);
            }
            reverse(all(SCCR));
            for (auto comp: SCCR){
                X.clear();
                Y.clear();
                for (auto v: comp)
                    for (auto e: G[v]){
                        int w = t[e];
                        if (sccId[v] == sccId[w])
                            X.pb(mp(mp(v, w), e));
                        if (sccId[v] != sccId[w])
                            Y.pb(mp(mp(v, w), e));
                    }
                SC[u].pb(mp(X, Y));
            }
            for (auto v: children[u]){
                vis[v] = 0;
                G[v].clear();
                G_inv[v].clear();
            }
            for (int i = 0; i < nl; i++){
                out[i] = 0;
                D[i].clear();
            }
        }
        reach.pop_back();
    }

    void decompose_and_sequence(){
        init();
        sequence.clear();
        vector<pair<PII,RegEx*> >edges;
        vector<RegEx*> P(m + 1), R(n + 1);
        vector<int> decompose_order = reach;
        auto cmp = [&] (int x, int y){
            return id[x] < id[y];
        };
        sort(all(decompose_order), cmp);
        for (auto u: decompose_order){
            if (children[u].empty())
                continue;
            for (auto v: children[u])
                for (auto e: nontree[v])
                    P[e] = eval_and_sequence(e);
            // edges.clear();
            // for (auto v: children[u])
            //     for (auto e: sibling[v])
            //         edges.pb(mp(mp(v, t[e]), P[e])); 
            // vector<PathSequence> Yu = eleminate(edges);
            vector<PathSequence> Yu;
            for (auto comp: SC[u]){
                edges.clear();
                for (auto edge: comp.ff)
                    edges.pb(mp(edge.ff, P[edge.ss]));
                for (auto path: eleminate(edges))
                    Yu.pb(path);
                for (auto edge: comp.ss)
                    Yu.pb(PathSequence(P[edge.ss], edge.ff.ff, edge.ff.ss));
            }
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
        for (int i = int(decompose_order.size()) - 2; i >= 0; i--){
            int v = decompose_order[i];
            sequence.pb(PathSequence(S[v], ancestor[v], v));
        }
    }
};


int main(){
    freopen("graph.txt", "r", stdin);
    int n, m;
    scanf("%d%d", &n, &m);
    vector<PII> edges;
    for (int i = 1; i <= m; i++){
        int u, v;
        scanf("%d%d", &u, &v);
        edges.pb(mp(u, v));
    }
    Tarjan T(n, m, edges);
    cout << *T.query(1, 4) << "\n";
}