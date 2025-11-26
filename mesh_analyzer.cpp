/*
mesh_analyzer.cpp

Extended C++17 single-file tool to parse MEDIT .mesh files and compute geometry descriptors.
Enhancements compared to the original version:
 - Parses permissively the common MEDIT sections: Vertices, Triangles, Quadrilaterals, Tetrahedra, Hexahedra.
 - Parses a flexible "Periodic" section (captures common patterns: pairs of vertex indices with
   translation vectors, or stores raw tokens if the exact format is unknown).
 - Computes per-element statistics for triangles, quads, tets and hexes:
     - element id (0-based), vertex indices, material tag
     - centroid, area (for surface elems) or volume (for volume elems)
     - edge lengths and simple quality metric where applicable
 - Optional JSON output (--json filename). If filename is "-" JSON is written to stdout.
 - Optional flag --per-element to include per-element arrays in JSON (can be large).
 - Optional flag --pretty to pretty-print JSON with indentation.
 - If --json is NOT provided, a human-readable summary is printed to stdout as before.

Build:
  g++ -std=c++17 -O2 mesh_analyzer.cpp -o mesh_analyzer

Run:
  ./mesh_analyzer [--json out.json] [--per-element] [--pretty] path/to/file.mesh

Notes:
 - Parser is permissive: indices are treated as 1-based and converted to 0-based.
 - Hexahedron volumes are approximated by a heuristic decomposition into tetrahedra.
 - Periodic parsing is best-effort: common format parsed (a b tx ty tz), otherwise tokens stored raw.
*/

#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x=0, y=0, z=0;
    Vec3()=default;
    Vec3(double a,double b,double c):x(a),y(b),z(c){}
    Vec3 operator+(const Vec3& o) const { return Vec3(x+o.x,y+o.y,z+o.z); }
    Vec3 operator-(const Vec3& o) const { return Vec3(x-o.x,y-o.y,z-o.z); }
    Vec3 operator*(double s) const { return Vec3(x*s,y*s,z*s); }
    Vec3 operator/(double s) const { return Vec3(x/s,y/s,z/s); }
};

static inline double dot(const Vec3& a, const Vec3& b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline Vec3 cross(const Vec3& a, const Vec3& b){
    return Vec3(a.y*b.z - a.z*b.y,
                a.z*b.x - a.x*b.z,
                a.x*b.y - a.y*b.x);
}
static inline double norm(const Vec3& a){ return sqrt(max(0.0, dot(a,a))); }

struct Triangle { array<int,3> v; int m=0; };
struct Quad { array<int,4> v; int m=0; };
struct Tetra { array<int,4> v; int m=0; };
struct Hex { array<int,8> v; int m=0; };

struct Edge {
    int a,b;
    Edge()=default;
    Edge(int x,int y){ if(x<y){a=x;b=y;} else {a=y;b=x;} }
    bool operator==(Edge const& o) const { return a==o.a && b==o.b; }
};
struct EdgeHash { size_t operator()(Edge const& e) const noexcept {
    return ((uint64_t)e.a<<32) ^ (uint64_t)e.b;
} };

// Periodic entry: try to capture (a,b,shift) if available, otherwise store raw tokens
struct PeriodicEntry {
    bool hasPair = false;
    int a=-1, b=-1;
    bool hasShift = false;
    Vec3 shift;
    vector<string> rawTokens;
};

// Simple 3x3 symmetric Jacobi eigen decomposition for PCA
void jacobiEigen3x3(double m[3][3], array<double,3>& evals, double evecs[3][3], int maxIter=50, double eps=1e-12){
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) evecs[i][j] = (i==j?1.0:0.0);
    double a[3][3];
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) a[i][j]=m[i][j];

    for(int iter=0; iter<maxIter; ++iter){
        int p=0,q=1;
        double maxOff = fabs(a[0][1]);
        if(fabs(a[0][2])>maxOff){ p=0;q=2; maxOff=fabs(a[0][2]); }
        if(fabs(a[1][2])>maxOff){ p=1;q=2; maxOff=fabs(a[1][2]); }
        if(maxOff < eps) break;
        double app = a[p][p];
        double aqq = a[q][q];
        double apq = a[p][q];
        double phi = 0.5 * atan2(2*apq, aqq - app);
        double c = cos(phi), s = sin(phi);
        for(int i=0;i<3;i++){
            double aip = a[i][p], aiq = a[i][q];
            a[i][p] = c*aip - s*aiq;
            a[p][i] = a[i][p];
            a[i][q] = s*aip + c*aiq;
            a[q][i] = a[i][q];
        }
        a[p][p] = c*c*app - 2*s*c*apq + s*s*aqq;
        a[q][q] = s*s*app + 2*s*c*apq + c*c*aqq;
        a[p][q] = a[q][p] = 0.0;
        for(int i=0;i<3;i++){
            double vip = evecs[i][p], viq = evecs[i][q];
            evecs[i][p] = c*vip - s*viq;
            evecs[i][q] = s*vip + c*viq;
        }
    }
    evals = {a[0][0], a[1][1], a[2][2]};
    array<int,3> idx = {0,1,2};
    sort(idx.begin(), idx.end(), [&](int i,int j){ return evals[i] > evals[j]; });
    array<double,3> evals_sorted;
    double evecs_sorted[3][3];
    for(int r=0;r<3;r++){
        evals_sorted[r] = evals[idx[r]];
        for(int c=0;c<3;c++) evecs_sorted[c][r] = evecs[c][idx[r]];
    }
    evals = evals_sorted;
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) evecs[i][j] = evecs_sorted[i][j];
}

double triangleArea(const Vec3& a, const Vec3& b, const Vec3& c){
    return 0.5 * norm(cross(b-a, c-a));
}

double tetraVolumeSigned(const Vec3& a,const Vec3& b,const Vec3& c,const Vec3& d){
    return dot(cross(b-a, c-a), d-a) / 6.0;
}

// Utility: join vector<string> with space
string joinTokens(const vector<string>& toks, size_t start=0){
    string s;
    for(size_t i=start;i<toks.size();++i){
        if(i!=start) s += " ";
        s += toks[i];
    }
    return s;
}

// JSON helpers (very small, manual)
string json_escape(const string& s){
    string out;
    out.reserve(s.size()+10);
    for(char c: s){
        switch(c){
            case '\"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if((unsigned char)c < 0x20){
                    char buf[8];
                    snprintf(buf, sizeof(buf), "\\u%04x", (unsigned char)c);
                    out += buf;
                } else out += c;
        }
    }
    return out;
}

string vec3_to_json(const Vec3& v, bool pretty=false, int indent=0){
    char buf[256];
    if(pretty){
        string pad(indent,' ');
        snprintf(buf, sizeof(buf), "{\n%s  \"x\": %.12g,\n%s  \"y\": %.12g,\n%s  \"z\": %.12g\n%s}", pad.c_str(), v.x, pad.c_str(), v.y, pad.c_str(), v.z, pad.c_str());
        return string(buf);
    } else {
        snprintf(buf, sizeof(buf), "{\"x\":%.12g,\"y\":%.12g,\"z\":%.12g}", v.x, v.y, v.z);
        return string(buf);
    }
}

int main(int argc, char** argv){
    if(argc < 2){
        cerr<<"Usage: "<<argv[0]<<" [--json out.json] [--per-element] [--pretty] path/to/file.mesh\n";
        return 1;
    }

    string filename;
    string json_out;
    bool per_element = false;
    bool pretty = false;
    vector<string> files;

    for(int i=1;i<argc;i++){
        string a = argv[i];
        if(a=="--json" && i+1<argc){ json_out = argv[++i]; }
        else if(a=="--per-element") per_element = true;
        else if(a=="--pretty") pretty = true;
        else files.push_back(a);
    }
    if(files.empty()){
        cerr<<"Error: missing .mesh file path\n";
        return 1;
    }
    filename = files[0];

    ifstream in(filename);
    if(!in){
        cerr<<"Failed to open "<<filename<<"\n";
        return 1;
    }

    vector<Vec3> vertices;
    vector<Triangle> triangles;
    vector<Quad> quads;
    vector<Tetra> tets;
    vector<Hex> hexes;
    vector<PeriodicEntry> periodic;

    // Read tokens, preserve comments removal
    string token;
    vector<string> toks;
    // Also preserve original line tokens to enable raw capture if needed
    // We'll do token-based parsing but allow capturing raw tokens per section.
    while(in >> token){
        if(token.size()>0 && (token[0]=='%' || token[0]=='#')){
            string rest; getline(in, rest);
            continue;
        }
        toks.push_back(token);
    }

    size_t p = 0;
    auto hasNext = [&](size_t k=0)->bool{ return p+k < toks.size(); };
    while(p < toks.size()){
        string s = toks[p++];
        if(s=="MeshVersionFormatted" || s=="MeshVersionFormatted:") { if(hasNext()) ++p; continue; }
        if(s=="Dimension") { if(hasNext()) ++p; continue; }

        if(s == "Vertices"){
            if(!hasNext()) break;
            int n = stoi(toks[p++]);
            vertices.reserve(n);
            for(int i=0;i<n && p+3<=toks.size(); ++i){
                double x = stod(toks[p++]);
                double y = stod(toks[p++]);
                double z = stod(toks[p++]);
                // optional reference integer: if next token is integer, consume it
                if(hasNext()){
                    // try parse as int; but since we only have tokens, accept and skip one token
                    // if it's clearly an integer (no dots) we'll skip it
                    string maybe = toks[p];
                    bool isInt = true;
                    for(char c: maybe) if(!(c=='-' || (c>='0'&&c<='9'))) { isInt=false; break; }
                    if(isInt) ++p;
                }
                vertices.emplace_back(x,y,z);
            }
            continue;
        }

        if(s == "Triangles"){
            if(!hasNext()) break;
            int n = stoi(toks[p++]);
            triangles.reserve(n);
            for(int i=0;i<n && p+4<=toks.size(); ++i){
                int a=stoi(toks[p++]) - 1;
                int b=stoi(toks[p++]) - 1;
                int c=stoi(toks[p++]) - 1;
                int ref = stoi(toks[p++]);
                triangles.push_back({{a,b,c}, ref});
            }
            continue;
        }

        if(s == "Quadrilaterals" || s == "Quadrangles" || s == "Quads"){
            if(!hasNext()) break;
            int n = stoi(toks[p++]);
            quads.reserve(n);
            for(int i=0;i<n && p+5<=toks.size(); ++i){
                int a=stoi(toks[p++]) - 1;
                int b=stoi(toks[p++]) - 1;
                int c=stoi(toks[p++]) - 1;
                int d=stoi(toks[p++]) - 1;
                int ref=stoi(toks[p++]);
                quads.push_back({{a,b,c,d}, ref});
            }
            continue;
        }

        if(s == "Tetrahedra" || s == "Tetra"){
            if(!hasNext()) break;
            int n = stoi(toks[p++]);
            tets.reserve(n);
            for(int i=0;i<n && p+5<=toks.size(); ++i){
                int a=stoi(toks[p++]) - 1;
                int b=stoi(toks[p++]) - 1;
                int c=stoi(toks[p++]) - 1;
                int d=stoi(toks[p++]) - 1;
                int ref=stoi(toks[p++]);
                tets.push_back({{a,b,c,d}, ref});
            }
            continue;
        }

        if(s == "Hexahedra" || s == "Hexahedron" || s == "Hex"){
            if(!hasNext()) break;
            int n = stoi(toks[p++]);
            hexes.reserve(n);
            for(int i=0;i<n && p+9<=toks.size(); ++i){
                array<int,8> v;
                for(int j=0;j<8;j++) v[j] = stoi(toks[p++]) - 1;
                int ref = stoi(toks[p++]);
                hexes.push_back({v, ref});
            }
            continue;
        }

        // Flexible Periodic parsing:
        // Common pattern: "Periodic" <n>
        // each entry: id_a id_b tx ty tz   (pair of vertex indices and translation vector)
        // But other files may have different content. We'll attempt to parse the common case,
        // otherwise capture raw tokens.
        if(s == "Periodic" || s == "PeriodicVertices" || s == "PeriodicEdges"){
            if(!hasNext()) break;
            int n = 0;
            bool haveCount = false;
            if(hasNext()){
                // Some files put an integer count
                string tok = toks[p];
                bool isInt = true;
                for(char c: tok) if(!(c=='-' || (c>='0'&&c<='9'))) { isInt=false; break; }
                if(isInt){
                    haveCount = true;
                    n = stoi(toks[p++]);
                }
            }
            if(!haveCount){
                // unknown count: try read until next known section by peeking ahead.
                // We'll try to read until we see a token that matches a known section label or EOF.
                size_t start = p;
                // find how many tokens until next known label or end
                size_t q = p;
                vector<string> known = {"Vertices","Triangles","Quadrilaterals","Tetrahedra","Hexahedra",
                                        "Edges","End","Boundary","Periodic","MeshVersionFormatted","Dimension"};
                while(q < toks.size() && find(known.begin(), known.end(), toks[q]) == known.end()) ++q;
                // consume all tokens as one raw block and store a single PeriodicEntry with rawTokens
                PeriodicEntry e;
                e.rawTokens.insert(e.rawTokens.end(), toks.begin()+p, toks.begin()+q);
                periodic.push_back(e);
                p = q;
                continue;
            } else {
                for(int i=0;i<n && p < toks.size(); ++i){
                    // try to parse: two ints followed by three doubles
                    PeriodicEntry e;
                    if(p < toks.size()){
                        // attempt parse two ints
                        if(p+1 < toks.size()){
                            bool okPair = true;
                            for(size_t k=0;k<2;k++){
                                string &tk = toks[p+k];
                                for(char c: tk) if(!(c=='-' || (c>='0'&&c<='9'))){ okPair=false; break; }
                            }
                            if(okPair){
                                e.a = stoi(toks[p++]) - 1;
                                e.b = stoi(toks[p++]) - 1;
                                e.hasPair = true;
                                // next try 3 doubles
                                if(p+2 < toks.size()){
                                    bool okShift = true;
                                    for(size_t k=0;k<3;k++){
                                        string &tk = toks[p+k];
                                        bool okNum = false;
                                        // try parse as double by checking presence of digits or dot or e/E
                                        for(char c: tk) if((c>='0'&&c<='9') || c=='.' || c=='e' || c=='E' || c=='+' || c=='-') { okNum = true; break; }
                                        if(!okNum) { okShift=false; break; }
                                    }
                                    if(okShift){
                                        double sx = stod(toks[p++]);
                                        double sy = stod(toks[p++]);
                                        double sz = stod(toks[p++]);
                                        e.shift = Vec3(sx,sy,sz);
                                        e.hasShift = true;
                                    }
                                }
                            } else {
                                // not ints: capture next token as raw
                                e.rawTokens.push_back(toks[p++]);
                                // also possibly capture more tokens until newline-like boundary; but tokenized input lacks lines.
                            }
                        } else {
                            // only one token remaining: capture as raw
                            e.rawTokens.push_back(toks[p++]);
                        }
                    }
                    periodic.push_back(e);
                }
            }
            continue;
        }

        // ignore other labels and advance conservatively if next token looks numeric count
        // This makes parser permissive: if label followed by an integer, skip that many tokens (best-effort).
        // But to be safe we just continue.
    }

    // Compute descriptors
    size_t nv = vertices.size(), ntri = triangles.size(), nquad = quads.size(), ntet = tets.size(), nhex = hexes.size();

    Vec3 minB(1e300,1e300,1e300), maxB(-1e300,-1e300,-1e300);
    Vec3 centroid(0,0,0);
    for(auto &v : vertices){
        minB.x = min(minB.x, v.x); minB.y = min(minB.y, v.y); minB.z = min(minB.z, v.z);
        maxB.x = max(maxB.x, v.x); maxB.y = max(maxB.y, v.y); maxB.z = max(maxB.z, v.z);
        centroid = centroid + v;
    }
    if(nv>0) centroid = centroid / (double)nv;

    unordered_set<Edge,EdgeHash> edges;
    auto add_edge = [&](int a,int b){ if(a<0||b<0||a>= (int)nv||b>= (int)nv) return; edges.insert(Edge(a,b)); };
    for(auto &t: triangles){
        add_edge(t.v[0],t.v[1]); add_edge(t.v[1],t.v[2]); add_edge(t.v[2],t.v[0]);
    }
    for(auto &q: quads){
        add_edge(q.v[0],q.v[1]); add_edge(q.v[1],q.v[2]); add_edge(q.v[2],q.v[3]); add_edge(q.v[3],q.v[0]);
    }
    for(auto &tet: tets){
        for(int i=0;i<4;i++) for(int j=i+1;j<4;j++) add_edge(tet.v[i], tet.v[j]);
    }
    for(auto &h: hexes){
        int idx[12][2] = {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
        for(int e=0;e<12;e++) add_edge(h.v[idx[e][0]], h.v[idx[e][1]]);
    }

    vector<double> edgeLengths;
    edgeLengths.reserve(edges.size());
    for(auto &e: edges){
        edgeLengths.push_back(norm(vertices[e.a] - vertices[e.b]));
    }
    double el_min=0, el_max=0, el_mean=0, el_std=0;
    if(!edgeLengths.empty()){
        el_min = *min_element(edgeLengths.begin(), edgeLengths.end());
        el_max = *max_element(edgeLengths.begin(), edgeLengths.end());
        double sum=0; for(auto d:edgeLengths) sum+=d;
        el_mean = sum/edgeLengths.size();
        double var=0; for(auto d:edgeLengths) var += (d-el_mean)*(d-el_mean);
        el_std = sqrt(var/edgeLengths.size());
    }

    double surfaceArea = 0;
    for(auto &t: triangles){
        surfaceArea += triangleArea(vertices[t.v[0]], vertices[t.v[1]], vertices[t.v[2]]);
    }
    for(auto &q: quads){
        surfaceArea += triangleArea(vertices[q.v[0]], vertices[q.v[1]], vertices[q.v[2]]);
        surfaceArea += triangleArea(vertices[q.v[0]], vertices[q.v[2]], vertices[q.v[3]]);
    }

    double volume = 0;
    for(auto &t: tets){
        double v = fabs(tetraVolumeSigned(vertices[t.v[0]], vertices[t.v[1]], vertices[t.v[2]], vertices[t.v[3]]));
        volume += v;
    }
    for(auto &h: hexes){
        if(nv==0) break;
        int v0 = h.v[0], v1=h.v[1], v2=h.v[2], v3=h.v[3], v4=h.v[4], v5=h.v[5], v6=h.v[6], v7=h.v[7];
        array<array<int,4>,5> tlist = {{
            {v0,v1,v3,v4},
            {v1,v2,v3,v6},
            {v1,v3,v4,v6},
            {v4,v5,v6,v1},
            {v1,v6,v7,v3}
        }};
        for(auto &tt : tlist){
            bool ok=true;
            for(int k=0;k<4;k++) if(tt[k]<0 || tt[k]>= (int)nv) ok=false;
            if(!ok) continue;
            double vv = fabs(tetraVolumeSigned(vertices[tt[0]], vertices[tt[1]], vertices[tt[2]], vertices[tt[3]]));
            volume += vv;
        }
    }

    // triangle qualities
    vector<double> triQual;
    triQual.reserve(triangles.size());
    for(auto &t: triangles){
        Vec3 a = vertices[t.v[0]];
        Vec3 b = vertices[t.v[1]];
        Vec3 c = vertices[t.v[2]];
        double A = triangleArea(a,b,c);
        double la = norm(b-a), lb = norm(c-b), lc = norm(a-c);
        double denom = la*la + lb*lb + lc*lc;
        double q = 0;
        if(denom>0) q = (4.0*sqrt(3.0)*A) / denom;
        triQual.push_back(q);
    }
    double tri_q_mean=0, tri_q_std=0, tri_q_min=0, tri_q_max=0;
    if(!triQual.empty()){
        tri_q_min = *min_element(triQual.begin(), triQual.end());
        tri_q_max = *max_element(triQual.begin(), triQual.end());
        double sum=0; for(auto v:triQual) sum+=v;
        tri_q_mean = sum/triQual.size();
        double var=0; for(auto v:triQual) var+=(v-tri_q_mean)*(v-tri_q_mean);
        tri_q_std = sqrt(var/triQual.size());
    }

    // PCA
    double cov[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    if(nv>0){
        for(auto &v: vertices){
            Vec3 d = v - centroid;
            cov[0][0] += d.x*d.x; cov[0][1] += d.x*d.y; cov[0][2] += d.x*d.z;
            cov[1][0] += d.y*d.x; cov[1][1] += d.y*d.y; cov[1][2] += d.y*d.z;
            cov[2][0] += d.z*d.x; cov[2][1] += d.z*d.y; cov[2][2] += d.z*d.z;
        }
        double inv = 1.0 / (double)nv;
        for(int i=0;i<3;i++) for(int j=0;j<3;j++) cov[i][j] *= inv;
    }
    array<double,3> eigvals = {0,0,0};
    double eigvecs[3][3];
    jacobiEigen3x3(cov, eigvals, eigvecs);

    double compactness = 0;
    if(surfaceArea>0) compactness = volume / pow(surfaceArea, 1.5);

    // Prepare per-element stats if requested
    struct ElemStat {
        string type;
        int id;
        vector<int> verts;
        int material;
        Vec3 centroid;
        double area_or_volume;
        vector<double> edge_lengths;
        double quality;
    };
    vector<ElemStat> elemStats;
    if(per_element){
        // Triangles
        for(size_t i=0;i<triangles.size();++i){
            auto &t = triangles[i];
            ElemStat es;
            es.type = "triangle";
            es.id = (int)i;
            es.verts = {t.v[0], t.v[1], t.v[2]};
            es.material = t.m;
            Vec3 a = vertices[t.v[0]], b = vertices[t.v[1]], c = vertices[t.v[2]];
            es.centroid = (a + b + c) / 3.0;
            es.area_or_volume = triangleArea(a,b,c);
            es.edge_lengths = { norm(b-a), norm(c-b), norm(a-c) };
            double denom = es.edge_lengths[0]*es.edge_lengths[0] + es.edge_lengths[1]*es.edge_lengths[1] + es.edge_lengths[2]*es.edge_lengths[2];
            es.quality = denom>0 ? (4.0*sqrt(3.0)*es.area_or_volume)/denom : 0.0;
            elemStats.push_back(move(es));
        }
        // Quads
        for(size_t i=0;i<quads.size();++i){
            auto &q = quads[i];
            ElemStat es;
            es.type = "quad";
            es.id = (int)i;
            es.verts = {q.v[0], q.v[1], q.v[2], q.v[3]};
            es.material = q.m;
            Vec3 a = vertices[q.v[0]], b = vertices[q.v[1]], c = vertices[q.v[2]], d = vertices[q.v[3]];
            es.centroid = (a + b + c + d) / 4.0;
            // area as two triangles
            double A1 = triangleArea(a,b,c);
            double A2 = triangleArea(a,c,d);
            es.area_or_volume = A1 + A2;
            es.edge_lengths = { norm(b-a), norm(c-b), norm(d-c), norm(a-d) };
            // quality: mean of two triangle qualities
            double denom1 = es.edge_lengths[0]*es.edge_lengths[0] + es.edge_lengths[1]*es.edge_lengths[1] + (norm(a-c)*norm(a-c));
            double q1 = denom1>0 ? (4.0*sqrt(3.0)*A1)/denom1 : 0.0;
            double denom2 = es.edge_lengths[2]*es.edge_lengths[2] + es.edge_lengths[3]*es.edge_lengths[3] + (norm(b-d)*norm(b-d));
            double q2 = denom2>0 ? (4.0*sqrt(3.0)*A2)/denom2 : 0.0;
            es.quality = 0.5*(q1+q2);
            elemStats.push_back(move(es));
        }
        // Tetras
        for(size_t i=0;i<tets.size();++i){
            auto &t = tets[i];
            ElemStat es;
            es.type = "tetra";
            es.id = (int)i;
            es.verts = {t.v[0], t.v[1], t.v[2], t.v[3]};
            es.material = t.m;
            Vec3 a = vertices[t.v[0]], b = vertices[t.v[1]], c = vertices[t.v[2]], d = vertices[t.v[3]];
            es.centroid = (a + b + c + d) / 4.0;
            es.area_or_volume = fabs(tetraVolumeSigned(a,b,c,d));
            // edge lengths
            vector<pair<int,int>> pairs = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
            for(auto pr: pairs){
                es.edge_lengths.push_back(norm(vertices[es.verts[pr.first]] - vertices[es.verts[pr.second]]));
            }
            // simple quality: scaled volume / sum(edge^2) (not a standard metric but gives relative indicator)
            double sumsq = 0; for(double L: es.edge_lengths) sumsq += L*L;
            if(sumsq>0) es.quality = pow(es.area_or_volume, 2.0/3.0) * 6.0 / sumsq; else es.quality = 0.0;
            elemStats.push_back(move(es));
        }
        // Hexahedra
        for(size_t i=0;i<hexes.size();++i){
            auto &h = hexes[i];
            ElemStat es;
            es.type = "hexa";
            es.id = (int)i;
            es.verts = {h.v[0],h.v[1],h.v[2],h.v[3],h.v[4],h.v[5],h.v[6],h.v[7]};
            es.material = h.m;
            Vec3 sum(0,0,0);
            bool ok=true;
            for(int vi=0;vi<8;vi++){
                if(h.v[vi] < 0 || h.v[vi] >= (int)nv) { ok=false; break; }
                sum = sum + vertices[h.v[vi]];
            }
            if(ok) es.centroid = sum / 8.0;
            else es.centroid = Vec3(0,0,0);
            // approximate volume by same decomposition used above
            double vol = 0;
            if(ok){
                int v0 = h.v[0], v1=h.v[1], v2=h.v[2], v3=h.v[3], v4=h.v[4], v5=h.v[5], v6=h.v[6], v7=h.v[7];
                array<array<int,4>,5> tlist = {{
                    {v0,v1,v3,v4},
                    {v1,v2,v3,v6},
                    {v1,v3,v4,v6},
                    {v4,v5,v6,v1},
                    {v1,v6,v7,v3}
                }};
                for(auto &tt : tlist){
                    vol += fabs(tetraVolumeSigned(vertices[tt[0]], vertices[tt[1]], vertices[tt[2]], vertices[tt[3]]));
                }
            }
            es.area_or_volume = vol;
            // edges
            int idxE[12][2] = {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
            for(int e=0;e<12;e++){
                int a = h.v[idxE[e][0]], b = h.v[idxE[e][1]];
                if(a>=0 && b>=0 && a<nv && b<nv) es.edge_lengths.push_back(norm(vertices[a]-vertices[b]));
            }
            // quality: heuristic based on volume vs edge lengths
            double sumsq=0; for(double L: es.edge_lengths) sumsq += L*L;
            if(sumsq>0) es.quality = pow(es.area_or_volume, 2.0/3.0) * 8.0 / sumsq; else es.quality = 0.0;
            elemStats.push_back(move(es));
        }
    }

    // Output: JSON if requested, else text report
    ostringstream report;
    report.setf(std::ios::fixed); report<<setprecision(6);

    report<<"Mesh Analysis Report for: "<<filename<<"\n\n";
    report<<"Counts:\n";
    report<<"  vertices:       "<<nv<<"\n";
    report<<"  triangles:      "<<ntri<<"\n";
    report<<"  quadrilaterals: "<<nquad<<"\n";
    report<<"  tetrahedra:     "<<ntet<<"\n";
    report<<"  hexahedra:      "<<nhex<<"\n\n";

    report<<"Bounding box:\n";
    report<<"  min: ("<<minB.x<<", "<<minB.y<<", "<<minB.z<<")\n";
    report<<"  max: ("<<maxB.x<<", "<<maxB.y<<", "<<maxB.z<<")\n";
    report<<"  extents: ("<<maxB.x-minB.x<<", "<<maxB.y-minB.y<<", "<<maxB.z-minB.z<<")\n\n";
    report<<"Centroid (vertex average): ("<<centroid.x<<", "<<centroid.y<<", "<<centroid.z<<")\n\n";
    report<<"Edge statistics:\n";
    report<<"  unique edges: "<<edges.size()<<"\n";
    if(!edgeLengths.empty()){
        report<<"  edge length min: "<<el_min<<"\n";
        report<<"  edge length max: "<<el_max<<"\n";
        report<<"  edge length mean: "<<el_mean<<"\n";
        report<<"  edge length std: "<<el_std<<"\n";
    } else {
        report<<"  (no edges found)\n";
    }
    report<<"\n";
    report<<"Surface area (triangles + quads triangulated): "<<surfaceArea<<"\n";
    report<<"Volume (tets exact; hex approximated): "<<volume<<"\n";
    report<<"Compactness V / (A^(3/2)) : "<<compactness<<"\n\n";
    report<<"Triangle quality (4*sqrt(3)*area / sum(edge^2))\n";
    if(!triQual.empty()){
        report<<"  mean: "<<tri_q_mean<<", std: "<<tri_q_std<<", min: "<<tri_q_min<<", max: "<<tri_q_max<<"\n";
    } else {
        report<<"  (no triangles)\n";
    }
    report<<"\n";
    report<<"PCA (covariance eigenvalues): ["<<eigvals[0]<<", "<<eigvals[1]<<", "<<eigvals[2]<<"]\n\n";

    // Periodic summary
    report<<"Periodic entries parsed: "<<periodic.size()<<"\n";
    for(size_t i=0;i<periodic.size();++i){
        auto &pe = periodic[i];
        report<<"  ["<<i<<"] ";
        if(pe.hasPair){
            report<<"pair: ("<<pe.a<<","<<pe.b<<")";
            if(pe.hasShift) report<<", shift: ("<<pe.shift.x<<","<<pe.shift.y<<","<<pe.shift.z<<")";
            report<<"\n";
        } else {
            if(!pe.rawTokens.empty()) report<<"raw: \""<<joinTokens(pe.rawTokens)<<"\"\n";
            else report<<"(empty)\n";
        }
    }

    // If JSON requested, build JSON structure
    if(!json_out.empty()){
        bool toStdout = (json_out == "-");
        ofstream fout;
        if(!toStdout){
            fout.open(json_out);
            if(!fout){
                cerr<<"Failed to open JSON output file: "<<json_out<<"\n";
                return 1;
            }
        }
        ostream &out = toStdout ? cout : fout;
        // Build JSON manually
        string indentStr = pretty ? "  " : "";
        auto writeIndent = [&](int lvl){
            if(pretty) for(int k=0;k<lvl;k++) out<<indentStr;
        };
        out<< (pretty ? "{\n" : "{");

        int lvl = 1;
        // counts
        writeIndent(lvl); out<<"\"counts\": "<<(pretty?"{\n":"{");
        if(pretty){ writeIndent(lvl+1); out<<"\"vertices\": "<<nv<<",\n"; writeIndent(lvl+1); out<<"\"triangles\": "<<ntri<<",\n"; writeIndent(lvl+1); out<<"\"quadrilaterals\": "<<nquad<<",\n"; writeIndent(lvl+1); out<<"\"tetrahedra\": "<<ntet<<",\n"; writeIndent(lvl+1); out<<"\"hexahedra\": "<<nhex<<"\n"; writeIndent(lvl); out<<"}"; }
        else { out<<"\"vertices\": "<<nv<<",\"triangles\": "<<ntri<<",\"quadrilaterals\": "<<nquad<<",\"tetrahedra\": "<<ntet<<",\"hexahedra\": "<<nhex<<"}"; }

        out<<(pretty?",\n":",");
        // bbox
        writeIndent(lvl); out<<"\"bbox\": "<<(pretty?("{\n"):"{");
        if(pretty){
            writeIndent(lvl+1); out<<"\"min\": "<<vec3_to_json(minB,true,lvl*2+2)<<",\n";
            writeIndent(lvl+1); out<<"\"max\": "<<vec3_to_json(maxB,true,lvl*2+2)<<"\n";
            writeIndent(lvl); out<<"}";
        } else {
            out<<"\"min\": "<<vec3_to_json(minB)<<",\"max\": "<<vec3_to_json(maxB)<<"";
            out<<"}";
        }
        out<<(pretty?",\n":",");

        // centroid
        writeIndent(lvl); out<<"\"centroid\": "<<(pretty?vec3_to_json(centroid,true,lvl*2):vec3_to_json(centroid));
        out<<(pretty?",\n":",");

        // edge stats
        writeIndent(lvl); out<<"\"edge_stats\": "<<(pretty?"{\n":"{");
        if(pretty){
            writeIndent(lvl+1); out<<"\"unique_edges\": "<<edges.size()<<",\n";
            writeIndent(lvl+1); out<<"\"min\": "<<el_min<<",\n";
            writeIndent(lvl+1); out<<"\"max\": "<<el_max<<",\n";
            writeIndent(lvl+1); out<<"\"mean\": "<<el_mean<<",\n";
            writeIndent(lvl+1); out<<"\"std\": "<<el_std<<"\n";
            writeIndent(lvl); out<<"}";
        } else {
            out<<"\"unique_edges\": "<<edges.size()<<",\"min\": "<<el_min<<",\"max\": "<<el_max<<",\"mean\": "<<el_mean<<",\"std\": "<<el_std<<"}";
        }
        out<<(pretty?",\n":",");

        // surface/volume/compactness
        writeIndent(lvl); out<<"\"surface_area\": "<<surfaceArea<<",\n";
        writeIndent(lvl); out<<"\"volume\": "<<volume<<",\n";
        writeIndent(lvl); out<<"\"compactness\": "<<compactness<<(pretty?",\n":",");

        // triangle quality summary
        writeIndent(lvl); out<<"\"triangle_quality\": "<<(pretty?"{\n":"{");
        if(pretty){
            writeIndent(lvl+1); out<<"\"mean\": "<<tri_q_mean<<",\n";
            writeIndent(lvl+1); out<<"\"std\": "<<tri_q_std<<",\n";
            writeIndent(lvl+1); out<<"\"min\": "<<tri_q_min<<",\n";
            writeIndent(lvl+1); out<<"\"max\": "<<tri_q_max<<"\n";
            writeIndent(lvl); out<<"}";
        } else {
            out<<"\"mean\": "<<tri_q_mean<<",\"std\": "<<tri_q_std<<",\"min\": "<<tri_q_min<<",\"max\": "<<tri_q_max<<"}";
        }
        out<<(pretty?",\n":",");

        // PCA
        writeIndent(lvl); out<<"\"pca\": "<<(pretty?"{\n":"{");
        if(pretty){
            writeIndent(lvl+1); out<<"\"eigenvalues\": ["<<eigvals[0]<<", "<<eigvals[1]<<", "<<eigvals[2]<<"],\n";
            writeIndent(lvl+1); out<<"\"eigenvectors\": [\n";
            for(int c=0;c<3;c++){
                writeIndent(lvl+2); out<<"["<<eigvecs[0][c]<<", "<<eigvecs[1][c]<<", "<<eigvecs[2][c]<<"]";
                if(c<2) out<<",\n"; else out<<"\n";
            }
            writeIndent(lvl+1); out<<"]\n";
            writeIndent(lvl); out<<"}";
        } else {
            out<<"\"eigenvalues\":["<<eigvals[0]<<","<<eigvals[1]<<","<<eigvals[2]<<"],\"eigenvectors\":[["<<eigvecs[0][0]<<","<<eigvecs[1][0]<<","<<eigvecs[2][0]<<"],["<<eigvecs[0][1]<<","<<eigvecs[1][1]<<","<<eigvecs[2][1]<<"],["<<eigvecs[0][2]<<","<<eigvecs[1][2]<<","<<eigvecs[2][2]<<"]]}";
        }
        out<<(pretty?",\n":",");

        // Periodic entries
        writeIndent(lvl); out<<"\"periodic\": "<<(pretty?"[\n":"[");
        if(pretty){
            for(size_t i=0;i<periodic.size();++i){
                auto &pe = periodic[i];
                writeIndent(lvl+1); out<<"{\n";
                if(pe.hasPair){
                    writeIndent(lvl+2); out<<"\"a\": "<<pe.a<<",\n";
                    writeIndent(lvl+2); out<<"\"b\": "<<pe.b<<",\n";
                    if(pe.hasShift){
                        writeIndent(lvl+2); out<<"\"shift\": "<<vec3_to_json(pe.shift,true,lvl*2+4)<<"\n";
                    } else {
                        writeIndent(lvl+2); out<<"\"shift\": null\n";
                    }
                } else {
                    writeIndent(lvl+2); out<<"\"raw\": \""<<json_escape(joinTokens(pe.rawTokens))<<"\"\n";
                }
                writeIndent(lvl+1); out<<"}";
                if(i+1<periodic.size()) out<<",\n"; else out<<"\n";
            }
            writeIndent(lvl); out<<"]";
        } else {
            for(size_t i=0;i<periodic.size();++i){
                auto &pe = periodic[i];
                out<<"{";
                if(pe.hasPair){
                    out<<"\"a\": "<<pe.a<<",\"b\": "<<pe.b<<",\"shift\": "<<vec3_to_json(pe.shift);
                } else {
                    out<<"\"raw\":\""<<json_escape(joinTokens(pe.rawTokens))<<"\"";
                }
                out<<"}";
                if(i+1<periodic.size()) out<<",";
            }
            out<<"]";
        }

        // Per-element details
        if(per_element){
            out<<(pretty?",\n":",");
            writeIndent(lvl); out<<"\"elements\": "<<(pretty?"[\n":"[");
            if(pretty){
                for(size_t i=0;i<elemStats.size();++i){
                    auto &es = elemStats[i];
                    writeIndent(lvl+1); out<<"{\n";
                    writeIndent(lvl+2); out<<"\"type\": \""<<es.type<<"\",\n";
                    writeIndent(lvl+2); out<<"\"id\": "<<es.id<<",\n";
                    writeIndent(lvl+2); out<<"\"material\": "<<es.material<<",\n";
                    writeIndent(lvl+2); out<<"\"verts\": [";
                    for(size_t k=0;k<es.verts.size();++k){ out<<es.verts[k]; if(k+1<es.verts.size()) out<<", "; }
                    out<<"],\n";
                    writeIndent(lvl+2); out<<"\"centroid\": "<<vec3_to_json(es.centroid,true,lvl*2+4)<<",\n";
                    writeIndent(lvl+2); out<<"\"area_or_volume\": "<<es.area_or_volume<<",\n";
                    writeIndent(lvl+2); out<<"\"edge_lengths\": [";
                    for(size_t k=0;k<es.edge_lengths.size();++k){ out<<es.edge_lengths[k]; if(k+1<es.edge_lengths.size()) out<<", "; }
                    out<<"],\n";
                    writeIndent(lvl+2); out<<"\"quality\": "<<es.quality<<"\n";
                    writeIndent(lvl+1); out<<"}";
                    if(i+1<elemStats.size()) out<<",\n"; else out<<"\n";
                }
                writeIndent(lvl); out<<"]";
            } else {
                for(size_t i=0;i<elemStats.size();++i){
                    auto &es = elemStats[i];
                    out<<"{\"type\":\""<<es.type<<"\",\"id\": "<<es.id<<",\"material\": "<<es.material<<",\"verts\":[";
                    for(size_t k=0;k<es.verts.size();++k){ out<<es.verts[k]; if(k+1<es.verts.size()) out<<","; }
                    out<<"],\"centroid\": "<<vec3_to_json(es.centroid)<<",\"area_or_volume\": "<<es.area_or_volume<<",\"edge_lengths\":[";
                    for(size_t k=0;k<es.edge_lengths.size();++k){ out<<es.edge_lengths[k]; if(k+1<es.edge_lengths.size()) out<<","; }
                    out<<"],\"quality\": "<<es.quality<<"}";
                    if(i+1<elemStats.size()) out<<",";
                }
                out<<"]";
            }
        }

        // finish JSON
        if(pretty) out<<"\n}\n"; else out<<"}\n";

        if(!toStdout) cout<<"Wrote JSON to "<<json_out<<"\n";
    } else {
        // print human-readable report
        cout<<report.str();
    }

    return 0;
}