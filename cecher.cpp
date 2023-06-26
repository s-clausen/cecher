/*
 Cecher: efficient computation of Čech persistence barcodes
 2023, Soenke Clausen
*/


#define EXPLICIT_CENTERS
//#define USE_RATIONALS
//#define USE_COEFFICIENTS
//#define INDICATE_PROGRESS
#define PRINT_PERSISTENCE_PAIRS
//#define DEBUG
//#define USE_ROBINHOOD_HASHMAP


#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <stack>

#include <CGAL/Exact_rational.h>
#include <CGAL/Interval_nt.h>

#ifdef USE_ROBINHOOD_HASHMAP

#include "robin-hood-hashing/src/include/robin_hood.h"

template <class Key, class T, class H, class E>
using hash_map = robin_hood::unordered_map<Key, T, H, E>;
template <class Key> using hash = robin_hood::hash<Key>;

#else
template <class Key, class T, class H, class E> using hash_map = std::unordered_map<Key, T, H, E>;
template <class Key> using hash = std::hash<Key>;
#endif

bool use_threshold(false);
int point_dim;

#ifdef USE_RATIONALS
typedef CGAL::Exact_rational ET;
typedef CGAL::Exact_rational PT;
#endif

#ifndef USE_RATIONALS
typedef CGAL::Exact_rational ET;
typedef CGAL::Interval_nt<> PT;
#endif

typedef std::vector<std::vector<ET>> mtx_ET;
typedef std::vector<std::vector<PT>> mtx_PT;

typedef int64_t index_t;
typedef uint16_t coefficient_t;

#ifdef INDICATE_PROGRESS
static const std::chrono::milliseconds time_step(40);
#endif

static const std::string clear_line("\r\033[K");
static const size_t num_coefficient_bits = 8;
static const index_t max_simplex_index =
        (index_t(1) << (8 * sizeof(index_t) - 1 - num_coefficient_bits)) - 1;

void check_overflow(index_t i) {
    if
#ifdef USE_COEFFICIENTS
        (i > max_simplex_index)
#else
            (i < 0)
#endif
        throw std::overflow_error("simplex index in filtration is larger than maximum index");
}

class table {
    std::vector<std::vector<index_t>> B;

public:
    table(index_t n, index_t k) : B(k + 1, std::vector<index_t>(n + 1, 0)) {
        for (index_t i = 0; i <= n; ++i) {
            B[0][i] = 1;
            for (index_t j = 1; j < std::min(i, k + 1); ++j)
                B[j][i] = B[j - 1][i - 1] + B[j][i - 1];
            if (i <= k) B[i][i] = 1;
            check_overflow(B[std::min(i >> 1, k)][i]);
        }
    }

    index_t operator()(index_t n, index_t k) const {
        return B[k][n];
    }
};

bool is_prime(const coefficient_t n) {
    if (!(n & 1) || n < 2) return n == 2;
    for (coefficient_t p = 3; p <= n / p; p += 2)
        if (!(n % p)) return false;
    return true;
}

std::vector<coefficient_t> multiplicative_inverse_vector(const coefficient_t m) {
    std::vector<coefficient_t> inverse(m);
    inverse[1] = 1;
    for (coefficient_t a = 2; a < m; ++a) inverse[a] = m - (inverse[m % a] * (m / a)) % m;
    return inverse;
}
index_t recomputation_count= 0;
index_t step_count= 0;
#ifdef DEBUG
auto start = std::chrono::steady_clock::now();
#endif
const char* filename = nullptr;
typedef std::pair<index_t,int> index_dim_t;

#ifdef USE_COEFFICIENTS
struct entry_t {
	index_t index : 8 * sizeof(index_t) - num_coefficient_bits;
	coefficient_t coefficient : num_coefficient_bits;
	entry_t(index_t _index, coefficient_t _coefficient)
		: index(_index), coefficient(_coefficient) {}
	entry_t(index_t _index) : index(_index), coefficient(0) {}
	entry_t() : index(0), coefficient(0) {}
};


entry_t make_entry(index_t i, coefficient_t c) { return entry_t(i, c); }
index_t get_index(const entry_t& e) { return e.index; }
index_t get_coefficient(const entry_t& e) { return e.coefficient; }
void set_coefficient(entry_t& e, const coefficient_t c) { e.coefficient = c; }

std::ostream& operator<<(std::ostream& stream, const entry_t& e) {
	stream << get_index(e) << ":" << get_coefficient(e);
	return stream;
}
#else

typedef index_t entry_t;
const index_t get_index(const entry_t& i) { return i; }
index_t get_coefficient(const entry_t& i) { return 1; }
entry_t make_entry(index_t _index, coefficient_t _value) { return entry_t(_index); }
void set_coefficient(entry_t& e, const coefficient_t c) {}

#endif

const entry_t& get_entry(const entry_t& e) { return e; }
typedef std::pair<PT, index_t> diameter_index_t;
PT get_diameter(const diameter_index_t& i) { return i.first; }
index_t get_index(const diameter_index_t& i) { return i.second; }
typedef std::pair<index_t, PT> index_diameter_t;
index_t get_index(const index_diameter_t& i) { return i.first; }
PT get_diameter(const index_diameter_t& i) { return i.second; }


template <typename T> struct ball_t : std::pair<T,index_dim_t>{
    ball_t(T _diameter, index_dim_t _support) : std::pair<T,index_dim_t> { _diameter,_support} {}
};

struct diameter_entry_t : std::tuple<PT, entry_t, index_dim_t> {
    using std::tuple<PT, entry_t, index_dim_t>::tuple;
    diameter_entry_t(const diameter_entry_t& p, coefficient_t _coefficient)
            : diameter_entry_t(get<0>(p), make_entry(get_index(get<1>(p)),_coefficient), get<2>(p)) {}
    diameter_entry_t(PT _diameter, index_t _index, coefficient_t _coefficient, index_dim_t _support_dim)
            : diameter_entry_t(_diameter, make_entry(_index, _coefficient), _support_dim) {}
    diameter_entry_t(PT _diameter, index_t _index, coefficient_t _coefficient)
            : diameter_entry_t(_diameter, make_entry(_index, _coefficient), index_dim_t(-1,-1)) {}
    diameter_entry_t(const diameter_index_t& p, coefficient_t _coefficient)
            : diameter_entry_t(get_diameter(p), make_entry(get_index(p),_coefficient), index_dim_t(get_index(p),1)) {}
    diameter_entry_t(const diameter_index_t& p)
            : diameter_entry_t(get_diameter(p), make_entry(get_index(p), 0),index_dim_t(get_index(p),1)) {}
    diameter_entry_t(const index_t& _index) : diameter_entry_t(0, _index, 0) {}
    diameter_entry_t(ball_t<PT> _ball, index_t _index, coefficient_t _coefficient)
            : diameter_entry_t(_ball.first, make_entry(_index, _coefficient),_ball.second) {}
#ifndef USE_RATIONALS
    diameter_entry_t(ball_t<ET> _ball, index_t _index, coefficient_t _coefficient)
            : diameter_entry_t(CGAL::to_interval(_ball.first), make_entry(_index, _coefficient), _ball.second) {}
#endif
};


const entry_t& get_entry(const diameter_entry_t& p) { return std::get<1>(p); }
entry_t& get_entry(diameter_entry_t& p) { return std::get<1>(p); }
const index_t get_index(const diameter_entry_t& p) { return get_index(get_entry(p)); }
const coefficient_t get_coefficient(const diameter_entry_t& p) { return get_coefficient(get_entry(p)); }
const PT& get_diameter(const diameter_entry_t& p) { return std::get<0>(p); }
const index_dim_t& get_support(const diameter_entry_t& p) { return std::get<2>(p); }
const index_t& get_supp_index(const diameter_entry_t& p) { return std::get<2>(p).first; }
const int& get_supp_dim(const diameter_entry_t& p) { return std::get<2>(p).second; }

void set_coefficient(diameter_entry_t& p, const coefficient_t c) {
    set_coefficient(get_entry(p), c);
}


template <typename T> struct upper_tri_mtx: std::vector<std::vector<T>>{};

template <typename T> std::vector<T> mtx_vector_product(const upper_tri_mtx<T>& M,const std::vector<T>& v){
    std::vector<T> result (v.size());
    for (int i=0; i<v.size(); i++) {
        result[i]=0;
        for (int j=0; j<v.size(); j++){
            if (i>=j) result[i]+=M[i][j] * v[j];
            else result[i]+=M[j][i] * v[j];
        }
    }
    return result;
}

template <typename T> upper_tri_mtx<T> outer_product(std::vector<T>& a,std::vector<T>& b){
    upper_tri_mtx<T> result;
    result.resize(a.size());
    for (int i=0;i< a.size(); i++) {
        for (int j=0;j< i+1; j++) {
            result[i].push_back(a[i]*b[j]);
        }
    }
    return result;
}

template <typename T> T inner_product(const std::vector<T>& a,const std::vector<T>& b){
    T result=0;
    for(int i=0; i<a.size(); i++) result+= a[i]*b[i];
    return result;
}

template <typename T> T compute_dist(const std::vector<T>& a,const std::vector<T>& b){
    std::vector<T> result(a.size());
    for (int i=0;i< a.size(); i++) result[i]=a[i]-b[i];
    return inner_product(result,result);
}

template <class Predicate>
index_t get_max(index_t top, const index_t bottom, const Predicate pred) {
    if (!pred(top)) {
        index_t count = top - bottom;
        while (count > 0) {
            index_t step = count >> 1, mid = top - step;
            if (!pred(mid)) {
                top = mid - 1;
                count -= step + 1;
            }
            else
                count = step;
        }
    }
    return top;
}


enum compressed_matrix_layout { LOWER_TRIANGULAR, UPPER_TRIANGULAR };

template <compressed_matrix_layout Layout> struct compressed_distance_matrix {
    std::vector<PT> distances;
    std::vector<PT*> rows;

    compressed_distance_matrix(std::vector<PT>&& _distances)
            : distances(std::move(_distances)), rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
        assert(distances.size() == size() * (size() - 1) / 2);
        init_rows();
    }

    template <typename DistanceMatrix>
    compressed_distance_matrix(const DistanceMatrix& mat)
            : distances(mat.size()* (mat.size() - 1) / 2), rows(mat.size()) {
        init_rows();

        for (size_t i = 1; i < size(); ++i)
            for (size_t j = 0; j < i; ++j) rows[i][j] = mat(i, j);
    }

    PT operator()(const index_t i, const index_t j) const;
    size_t size() const { return rows.size(); }
    void init_rows();
};

typedef compressed_distance_matrix<LOWER_TRIANGULAR> distance_matrix;

template <> void distance_matrix::init_rows() {
    PT* pointer = &distances[0];
    for (size_t i = 1; i < size(); ++i) {
        rows[i] = pointer;
        pointer += i;
    }
}

template <>
PT distance_matrix::operator()(const index_t i, const index_t j) const {
    return i == j ? 0 : i < j ? rows[j][i] : rows[i][j];
}

#ifdef INDICATE_PROGRESS
std::chrono::steady_clock::time_point next = std::chrono::steady_clock::now() + time_step;
#endif

struct euclidean_distance_matrix {

    mtx_PT points_PT;

    euclidean_distance_matrix(mtx_PT&& _points) : points_PT(std::move(_points)){
        for (auto p : points_PT) { assert(p.size() == points.front().size()); }
    }

    PT operator()(const index_t i, const index_t j) const {
        assert(i < points.size());
        assert(j < points.size());
#ifdef INDICATE_PROGRESS
        if (std::chrono::steady_clock::now() > next) {
            std::cerr << clear_line << "constructing distance matrix (processing "
                      << i << "/" << points_PT.size() << " vertices)" << std::flush;
            next = std::chrono::steady_clock::now() + time_step;
        }
#endif
        PT distance=compute_dist<PT>(points_PT[i],points_PT[j]);
        if(distance<=0) throw std::invalid_argument("algorithm cannot handle sets with vanishing distances");

        return distance;
    }

    size_t size() const { return points_PT.size(); }
};

class union_find {
    std::vector<index_t> parent;
    std::vector<uint8_t> rank;

public:
    union_find(const index_t n) : parent(n), rank(n, 0) {
        for (index_t i = 0; i < n; ++i) parent[i] = i;
    }

    index_t find(index_t x) {
        index_t y = x, z;
        while ((z = parent[y]) != y) y = z;
        while ((z = parent[x]) != y) {
            parent[x] = y;
            x = z;
        }
        return z;
    }

    void link(index_t x, index_t y) {
        if ((x = find(x)) == (y = find(y))) return;
        if (rank[x] > rank[y])
            parent[y] = x;
        else {
            parent[x] = y;
            if (rank[x] == rank[y]) ++rank[y];
        }
    }
};

template <typename T> T begin(std::pair<T, T>& p) { return p.first; }
template <typename T> T end(std::pair<T, T>& p) { return p.second; }

template <typename ValueType> class compressed_sparse_matrix {
    std::vector<size_t> bounds;
    std::vector<ValueType> entries;

    typedef typename std::vector<ValueType>::iterator iterator;
    typedef std::pair<iterator, iterator> iterator_pair;

public:
    size_t size() const { return bounds.size(); }

    iterator_pair subrange(const index_t index) {
        return { entries.begin() + (index == 0 ? 0 : bounds[index - 1]),
                 entries.begin() + bounds[index] };
    }

    void append_column() { bounds.push_back(entries.size()); }

    void push_back(const ValueType e) {
        assert(0 < size());
        entries.push_back(e);
        ++bounds.back();
    }
};



template <typename T> class min_ball {

    const std::vector<std::vector<T>>& points;
    const table& binomial_coeff;
    const distance_matrix& dist;

    bool is_assembly,init_coboundary;
    int idx, dim, m;
    std::vector<index_t> vertices;
    std::vector<int> indices;

    std::vector<std::vector<upper_tri_mtx<T>>> cache;
    upper_tri_mtx<T> prev_CM_inv, CM_inv;
    T diameter, radius;
    std::vector<T> center;

    upper_tri_mtx<T> simplex_CM_inv;
    std::vector<int> simplex_indices;
    std::vector<T> simplex_center;
    T simplex_diameter, simplex_radius;

    std::vector<index_dim_t> known_intvs;
    bool new_predicates;
    std::vector<T> radii;
    std::vector<std::vector<T>> centers;

    T half;

public:

    min_ball(const std::vector<std::vector<T>>& _points,
             const table& _binomial_coeff, const distance_matrix& _dist) :
            points(_points), binomial_coeff(_binomial_coeff), dist(_dist)
    {half = 1; half= half/2; center.resize(point_dim);}

    T get_dist(const index_t& a, const index_t& b);

    void throw_to_exact_type();

    std::vector<int> colexicographic_subsequent(const int max) {
        int size = indices.size();
        for (int l=0; l<size; l++){

            if (l+1< size && indices[l]<indices[l+1]-1) {
                indices[l]+=1;
                for (int i=0; i<l; i++) indices[i]=i;
                return indices;
            }
            else if (l+1==size && indices[l]<max-1) {
                for (int i=0; i<size-1; i++) indices[i]=i;
                indices[l]+=1;
                for (int i=0; i<l; i++) indices[i]=i;
                return indices;
            }
        }
        for (int l=0; l<size; l++) indices[l]=l;
        indices.push_back(size);
        return indices;
    }


    index_t compute_support_index() {
        index_t support=0;
        std::vector<int> order(dim+1);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),[&](int a, int b) {return vertices[indices[a]]<vertices[indices[b]];});
        for (int i=0; i<dim+1; i++) support+=binomial_coeff(vertices[indices[order[i]]],i+1);
        return support;
    }


    T compute_radius(const upper_tri_mtx<T>& p) {return -p[0][0]/2;}


    std::vector<T> compute_center(const upper_tri_mtx<T>& CM_inv, const int& max) {
        std::vector<T> result (point_dim,0);
        for (int i = 0; i < point_dim; i++) {
            for (int k = 0; k < max; k++)
                result[i] += points[vertices[indices[k]]][i] * CM_inv[k+1][0];
        }
        return result;
    }


    void set(const bool& _is_assembly,const bool& _init_coboundary) {
        is_assembly=_is_assembly;
        init_coboundary=_init_coboundary;
    }


    void reset(const std::vector<index_t>& _vertices, const int& j, const bool simplex_cache=false) {
        vertices=_vertices;
        if (j!=-1) vertices.push_back(j);
        known_intvs.resize(0);

        if (!simplex_cache) {
            m=vertices.size();
            cache.clear();
            cache.resize(m);
            idx=0;
            dim=1;
            indices.resize(2);
            indices[0]=0;
            indices[1]=1;
        }
    }


    void update_indices_idx_dim(bool pushing_a_vertex=false) {
        indices = colexicographic_subsequent(m);
        if (idx<binomial_coeff(m,dim+1)-1) idx++;
        else {
            dim++;
            if (pushing_a_vertex) {
                idx=binomial_coeff(m-1,dim+1);
                indices[indices.size()-1]= m-1;
            }
            else idx=0;
        }
    }


    upper_tri_mtx<T> init_CM_inv(){
        upper_tri_mtx<T> result;
        result.resize(3);
        result[0].push_back(-get_dist(vertices[indices[0]],vertices[indices[1]])/2);
        for (int i=1;i<3; i++) result[i].push_back(half);
        T auxiliary = 1;
        auxiliary=auxiliary/(2*get_dist(vertices[indices[0]],vertices[indices[1]]));
        result[2].push_back(auxiliary);
        for (int i=1;i<3; i++) result[i].push_back(-auxiliary);
        return result;
    }


    upper_tri_mtx<T> compute_CM_inv(const upper_tri_mtx<T>& prev_CM_inv,const std::vector<T>& CM_col) {
        upper_tri_mtx<T> result;

        std::vector<T> my = mtx_vector_product(prev_CM_inv, CM_col);
        T z = -1 * inner_product(CM_col, my);
        if (z==0) return result;
        T z_inv = 1/z;

        result = prev_CM_inv;
        std::vector<T> CM_inv_col = my;
        for (int i=0; i<CM_inv_col.size(); i++) CM_inv_col[i] *= -z_inv;
        upper_tri_mtx<T> product = outer_product(CM_inv_col, my);
        for (int i=0;i< result.size(); i++) {
            for (int j=0;j< i+1; j++) result[i][j]-=product[i][j];
        }
        CM_inv_col.push_back(z_inv);
        result.push_back(CM_inv_col);

        return result;
    }


    std::vector<T> next_CM_column(const index_t& vertex, const int& max) {
        std::vector<T> CM_col{1};
        for (int i = 0; i < max; i++) CM_col.push_back(get_dist(vertices[indices[i]], vertex));
        return CM_col;
    }


    void next_circumsphere() {
        if (dim==1) diameter=get_dist(vertices[indices[0]],vertices[indices[1]]);
        else {
            prev_CM_inv= cache[dim-2][idx- binomial_coeff(indices[dim],dim+1)];

            if (prev_CM_inv.size() == 0) {
                CM_inv.resize(0);
                return;
            }

            std::vector<T> CM_col=next_CM_column(vertices[indices[dim]],dim);

            CM_inv=compute_CM_inv(prev_CM_inv,CM_col);
            if(CM_inv.size()>0) radius=compute_radius(CM_inv);
        }
    }


    void get_predicates() {
        new_predicates=false;

        centers.clear();
        radii.clear();

        for (int d=0; d<dim-1;d++){

            int facet_idx=idx- binomial_coeff(indices[dim],dim+1);
            upper_tri_mtx<T> CM_inv = cache[dim-d-2][facet_idx];
            std::vector<T> center = compute_center(CM_inv,dim-d);

            if (compute_dist(points[vertices[indices[dim-d]]], center)
                == compute_radius(CM_inv)) {

                std::vector<T> push_point (point_dim);
                for (int l=0; l< point_dim; l++) {
                    push_point[l]= points[vertices[indices[dim-d]]][l]/2
                                   + points[vertices[indices[dim-d-1]]][l]/2;
                }
                std::vector<T> CM_col {1};
                for (int i=0; i<dim-d; i++) {
                    CM_col.push_back(compute_dist(points[vertices[indices[i]]],push_point));
                }
                CM_inv=compute_CM_inv(CM_inv,CM_col);

                center = compute_center(CM_inv,dim-d);
                for (int i = 0; i < point_dim; i++) center[i] += push_point[i] * CM_inv[dim-d+1][0];
            }
            centers.push_back(center);
            radii.push_back(compute_radius(CM_inv));

            if (d<dim-2) facet_idx-=binomial_coeff(indices[dim-d-1],dim-d);
        }
    }


    bool is_in_sphere(const index_t& vertex) {
        if (dim==1) {
            T orientation= -diameter;
            for (int i=0; i<2; i++) orientation+=get_dist(vertices[indices[i]],vertex);
            return (orientation<0 || orientation==0 && (vertex>vertices[indices[0]] || vertex>vertices[indices[1]]));
        }
#ifndef EXPLICIT_CENTERS
        std::vector<T> CM_col=next_CM_column(vertex,dim+1);
        std::vector<T> my=mtx_vector_product(CM_inv,CM_col);
        T minus_z=inner_product(CM_col,my);
        if (my[0]!=0 && minus_z!=0) return (my[0]/minus_z < 0);

        center=compute_center(CM_inv,dim+1);
#endif
        T margin = compute_dist(points[vertex], center);
        if (margin > radius) return false;
        else if (margin < radius) return true;

#ifndef USE_RATIONALS
        throw_to_exact_type();
#endif
        if (new_predicates) get_predicates();

        int dim=indices.size()-1;
        for (int d=0; d<dim;d++)  {

            if(d==dim-1) return (vertex>vertices[indices[dim-d]]);

            if (vertex<vertices[indices[dim-d]]) return false;

            margin=compute_dist(points[vertex], centers[d]);
            if (margin<radii[d]) return false;
            if (margin>radii[d]) return true;
        }
    }


    bool is_optimal() {
        for (int k=0; k< indices.size();k++) if (CM_inv[k+1][0]<0) return false;
        return true;
    }


    bool is_valid(){
        new_predicates=true;
#ifdef EXPLICIT_CENTERS
        if (dim>1) center=compute_center(CM_inv,dim+1);
#endif
        std::vector<int> in_sphere;
        int k=0;
        for(int i=0; i<vertices.size(); i++){
            if (k<=indices.back() && i==indices[k]) k++;
            else if (is_in_sphere(vertices[i])) {
                in_sphere.push_back(i);
            }
        }

        if (in_sphere.size()==vertices.size()-dim-1) return true;

        std::vector<std::vector<int>> intv;
        intv.push_back(indices);
        for (int i=0; i<in_sphere.size(); i++) {
            for (int j=0;2*j<intv.size();j++) {

                intv.push_back(intv[j]);
                intv.back().push_back(in_sphere[i]);
                std::sort(intv.back().begin(),intv.back().end());
                int idx=0;
                for (int k=0; k<intv.back().size(); k++) idx+=binomial_coeff(intv.back()[k],k+1);
                known_intvs.push_back(index_dim_t {idx,intv.back().size()-1});
            }
        }
        return false;
    }


    bool has_zero_cofacet(){
        new_predicates=true;
#ifdef EXPLICIT_CENTERS
        if (dim>1) center=compute_center(CM_inv,dim+1);
#endif
        std::vector<index_t> tmp = vertices;
        int k=0;
        std::sort(tmp.begin(), tmp.end());
        for(index_t i=0; i<points.size(); i++){
            if (k<=tmp.back() && i==tmp[k]) k++;
            else if (is_in_sphere(i)) return true;
        }
        return false;
    }


    void store_circumsphere(const bool all_spheres = true){
        if (dim==1 && (all_spheres || indices[dim]<m-1)) cache[dim-1].push_back(init_CM_inv());
        if (dim>1 && (all_spheres || indices[dim]<m-1)) cache[dim-1].push_back(CM_inv);
    }


    ball_t<T> construct(const std::vector<index_t>& _vertices, const index_t j=-1) {
        reset(_vertices, j);
        std::reverse(vertices.begin(),vertices.end());
        std::sort(vertices.begin(),vertices.end(),std::greater<index_t>());

        while (true){
            next_circumsphere();

            if(dim==m-1 || (dim==1 && is_valid()) || (dim>1 && CM_inv.size()>0 && is_optimal() &&
                      std::find(known_intvs.begin(), known_intvs.end(), index_dim_t{idx,dim})==known_intvs.end() && is_valid())){

                index_t support=0;
                if (is_assembly && (dim < m-1 || has_zero_cofacet())) support=-1;
                else support=compute_support_index();

                return ball_t<T>((dim>1)? radius*4: diameter,index_dim_t {support, dim});
            }

            store_circumsphere(false);
            update_indices_idx_dim();
        }
    }


    void construct_simplex_cache(const std::vector<index_t>& _vertices) {
        reset(_vertices, -1);
        bool ball_found =false;

        while (true){
            next_circumsphere();

            if(!ball_found && (dim==m-1 || (dim==1 && is_valid()) || (dim>1 && CM_inv.size()>0 && is_optimal() &&
                   std::find(known_intvs.begin(), known_intvs.end(), index_dim_t{idx,dim})==known_intvs.end() && is_valid()))){

#ifdef EXPLICIT_CENTERS
                if (dim>1) center=compute_center(CM_inv,dim+1);
#endif
                simplex_CM_inv = CM_inv;
                simplex_indices = indices;
                simplex_center = center;
                simplex_diameter= diameter;
                simplex_radius= radius;
                ball_found=true;
            }

            store_circumsphere();
            if (dim==m-1) {
                m++;
                return;
            }
            update_indices_idx_dim();
        }
    }


    ball_t<T> construct_cofacet(const std::vector<index_t>& _vertices, const index_t j=-1) {

        reset(_vertices, j, true);

        indices=simplex_indices;
        dim=indices.size()-1;
        idx=0;
        diameter=simplex_diameter;
        center=simplex_center;
        radius=simplex_radius;
#ifndef EXPLICIT_CENTERS
        CM_inv=simplex_CM_inv;
#endif
        if (!init_coboundary && is_in_sphere(j)){
            index_t support=0;
            if (is_assembly) support=-1;
            else support=compute_support_index();
            return ball_t<T>((dim>1)? radius*4: diameter, index_dim_t {support, dim});
        }
        indices.resize(2);
        indices[0]=0;
        indices[1]=m-1;
        dim=1;
        idx=binomial_coeff(m-1,2);

        while (true){
            if(dim==1 ||(dim>1 && std::find(known_intvs.begin(), known_intvs.end(), index_dim_t{idx,dim})==known_intvs.end())){
                next_circumsphere();

		if(dim==m-1 || (dim==1 && is_valid()) || (dim>1 && CM_inv.size()>0 && is_optimal() && is_valid())){

                    index_t support=0;
                    if (is_assembly && (dim < m-1 || has_zero_cofacet())) support=-1;
                    else support=compute_support_index();

                    return ball_t<T>((dim>1)? radius*4: diameter,index_dim_t {support, dim});
                }
            }

            update_indices_idx_dim(true);
        }
    }


    upper_tri_mtx<T> circumsphere(const std::vector<index_t>& _vertices) {
        reset(_vertices, -1);

        while (true){
            next_circumsphere();
            if(dim==m-1) return (dim==1)? init_CM_inv() : CM_inv;
            store_circumsphere();
            dim++;
            indices.push_back(dim);
        }
    }
};


template<> PT min_ball<PT>::get_dist(const index_t& a, const index_t& b){
    return dist(a,b);
}
#ifndef USE_RATIONALS
template<> ET min_ball<ET>::get_dist(const index_t& a, const index_t& b){
    return compute_dist(points[a],points[b]);
}
#endif

#ifndef USE_RATIONALS
template<> void min_ball<CGAL::Interval_nt<>>::throw_to_exact_type(){
    throw CGAL::Uncertain_conversion_exception("");
}
#endif
template<> void min_ball<CGAL::Exact_rational>::throw_to_exact_type(){}

template <typename T> class min_2_ball {

    const table &binomial_coeff;
    const distance_matrix &dist;
    std::vector<index_t> vertices;
    diameter_entry_t simplex;
    bool init_coboundary;

public:

    min_2_ball(const table &_binomial_coeff, const distance_matrix &_dist) :
            binomial_coeff(_binomial_coeff), dist(_dist) {}

    void set(const diameter_entry_t& _simplex,const bool& _init_coboundary) {
        init_coboundary=_init_coboundary;
        simplex=_simplex;
    }

    index_t compute_1_support(const index_t& a, const index_t& b) {
        if (a < b) return a+binomial_coeff(b,2);
        else return b+binomial_coeff(a,2);
    }

    bool is_smallest_index(const index_t& index) {
        for (int i=0;i<vertices.size();i++) if (vertices[i]<index) return false;
        return true;
    }

    ball_t<PT> construct(const std::vector<index_t>& _vertices, const index_t& j,const index_t& cofacet_index){
        vertices=_vertices;
        vertices.push_back(j);
        index_t support;

        std::vector<PT> edges {dist(j,vertices[1]),dist(j,vertices[0]),get_diameter(simplex)};

        PT cofacet_diameter = edges[2-init_coboundary];
        int max_edge_index = 2-init_coboundary;

        for (int i = 0; i < 2-init_coboundary; i++) {
            if (edges[i] > cofacet_diameter) {
                cofacet_diameter = edges[i];
                max_edge_index = i;
            }
        }

        PT orientation = -edges[max_edge_index];
        for (int i = 0; i < 3; i++) if (i != max_edge_index) orientation += edges[i];

        if (orientation == 0 && is_smallest_index(vertices[max_edge_index])) max_edge_index=-1;
        else if (orientation > 0) {
            PT auxiliary = 4;
            for (int i = 0; i < 3; i++) if (i != max_edge_index) auxiliary *= edges[i];
            cofacet_diameter = (edges[max_edge_index]*auxiliary) / (auxiliary-(orientation*orientation));
            max_edge_index=-1;
        }

        if (max_edge_index==-1) support= cofacet_index;
        else if (max_edge_index==2) support= get_index(simplex);
        else support = compute_1_support(j,vertices[(max_edge_index==0)? 1 : 0]);

        return ball_t<PT>(cofacet_diameter,index_dim_t {support, (max_edge_index==-1)? 2 : 1});
    }
};



template <typename DistanceMatrix> class cecher {

    const mtx_PT points_PT;
    const mtx_ET points_ET;
    const DistanceMatrix dist;
    const index_t n, dim_max;
    const PT threshold;
    const PT ratio;
    const coefficient_t modulus;
    const table binomial_coeff;
    const std::vector<coefficient_t> multiplicative_inverse;
    mutable std::vector<diameter_entry_t> cofacet_entries;
    mutable std::vector<index_t> vertices;

    min_ball<ET> min_ball_ET;

    struct entry_hash {
        std::size_t operator()(const entry_t& e) const { return hash<index_t>()(::get_index(e)); }
    };

    struct equal_index {
        bool operator()(const entry_t& e, const entry_t& f) const {
            return ::get_index(e) == ::get_index(f);
        }
    };

    typedef hash_map<entry_t, size_t, entry_hash, equal_index> entry_hash_map;

public:
    cecher(DistanceMatrix&& _dist, index_t _dim_max, PT _threshold, PT _ratio,
           coefficient_t _modulus, mtx_PT&& _points_PT, mtx_ET&& _points_ET)
            : dist(std::move(_dist)), n(dist.size()),
              dim_max(std::min(_dim_max, index_t(dist.size()-2))), threshold(_threshold),
              ratio(_ratio), modulus(_modulus), binomial_coeff(n, dim_max+2),
              multiplicative_inverse(multiplicative_inverse_vector(_modulus)),
              points_PT(std::move(_points_PT)), points_ET(std::move(_points_ET)),
              min_ball_ET(min_ball<ET>(points_ET,binomial_coeff,dist)) {}


    void print_pairs(bool upper_finite, PT lower, PT upper) {
        bool print=false;

        try {if (lower<upper) print=true;}
        catch(CGAL::Uncertain_conversion_exception){
            print=true;
        }
        if (!upper_finite) std::cout<<" ["<< std::sqrt(CGAL::to_double(lower))/2<<",)"<< std::endl;
        else std::cout<<" ["<<std::sqrt(CGAL::to_double(lower))/2<<","<<std::sqrt(CGAL::to_double(upper))/2<<")"<< std::endl;
    }


    static index_t get_max_vertex(const index_t idx, const index_t k, const index_t n, const table& binomial_coeff) {
        return get_max(n, k - 1, [&](index_t w) -> bool { return (binomial_coeff(w, k) <= idx); });
    }

    template <typename OutputIterator> static
    OutputIterator get_simplex_vertices(index_t idx, const index_t dim, index_t n,
                                        OutputIterator out, const table& binomial_coeff) {
        --n;
        for (index_t k = dim + 1; k > 1; --k) {
            n = get_max_vertex(idx, k, n, binomial_coeff);
            *out++ = n;
            idx -= binomial_coeff(n, k);
        }
        *out = idx;
        return out;
    }

    class simplex_coboundary_enumerator;


    diameter_entry_t get_zero_apparent_facet(const diameter_entry_t simplex, const index_t dim) {
        if (get_supp_dim(simplex)<dim) {

            if (get_supp_dim(simplex) ==dim-1) return diameter_entry_t(get_diameter(simplex), get_supp_index(simplex), 1);
            else {

                std::vector<index_t> support(get_supp_dim(simplex)+1), vertices(dim+1);
                get_simplex_vertices(get_index(simplex), dim, n, vertices.rbegin(), binomial_coeff);
                get_simplex_vertices(get_supp_index(simplex), get_supp_dim(simplex), n, support.rbegin(), binomial_coeff);

                index_t facet_index=0, highest_vertex=0, k=0;
                for (int i=0;i<vertices.size();i++) {
                    if (k<support.size() && vertices[i]==support[k]) k++;
                    else highest_vertex=vertices[i];
                }

                std::vector<index_t> facet;
                for (int i=0;i<dim+1;i++) if (vertices[i]!=highest_vertex) facet.push_back(vertices[i]);
                std::sort(facet.begin(), facet.end());
                for (int i=0; i<facet.size(); i++) facet_index+=binomial_coeff(facet[i],i+1);

                return diameter_entry_t(get_diameter(simplex), facet_index, 1);
            }
        }
        return diameter_entry_t(-1);
    }


    template <typename T> bool is_lower(const T& diameter,const T& threshold){
        try{return (diameter<=threshold);}
        catch(CGAL::Uncertain_conversion_exception) {return true;}
    }


    void assemble_columns_to_reduce(std::vector<diameter_index_t>& simplices,
                                    std::vector<diameter_entry_t>& columns_to_reduce,
                                    entry_hash_map& pivot_column_index, index_t dim) {
#ifdef DEBUG
        auto end = std::chrono::steady_clock::now();
        double elapsed_seconds = std::chrono::duration_cast< std::chrono::duration<double> >(end - start).count();
        std::cout << std::endl;
        std::cout << "### " << dim-1 << "-persistence in: " << elapsed_seconds << "s ###" <<  '\n';
        std::cout << "recomputations: " << recomputation_count <<  '\n';
        std::cout << "steps: " << step_count <<  '\n';
#endif

#ifdef INDICATE_PROGRESS
        std::cerr << clear_line << "assembling columns" << std::flush;
        std::chrono::steady_clock::time_point next = std::chrono::steady_clock::now() + time_step;
#endif
        columns_to_reduce.clear();
        std::vector<diameter_index_t> next_simplices;

        static simplex_coboundary_enumerator cofacets(*this);
        for (diameter_index_t& simplex : simplices) {
            cofacets.set_simplex(diameter_entry_t(simplex, 1), dim - 1, true);

            while (cofacets.has_next(false)) {

#ifdef INDICATE_PROGRESS
                if (std::chrono::steady_clock::now() > next) {
                    std::cerr << clear_line << "assembling"
                              << " columns (processing " << std::distance(&simplices[0], &simplex)
                              << "/" << simplices.size() << " simplices)" << std::flush;
                    next = std::chrono::steady_clock::now() + time_step;
                }
#endif

                auto cofacet = cofacets.next(true, false);

                if (!use_threshold || is_lower(get_diameter(cofacet),threshold)) {
                    if (dim < dim_max) next_simplices.push_back({ get_diameter(cofacet), get_index(cofacet) });
                    if (get_supp_index(cofacet)!=-1&&pivot_column_index.find(get_entry(cofacet))==pivot_column_index.end()){
                        columns_to_reduce.push_back(cofacet);
                    }
                }
            }
        }

        if (dim < dim_max) simplices.swap(next_simplices);

#ifdef DEBUG
#ifdef INDICATE_PROGRESS
        std::cerr << clear_line << std::flush;
#endif
        std::cout << "assembled " <<columns_to_reduce.size() << "/"
                  << binomial_coeff(n,dim+1) <<" "<< dim << "-simplices" <<std::endl;
        std::cout << std::endl;
#endif

        auto cmp = [&]
                (const diameter_entry_t& a, const diameter_entry_t& b) {
            return cofacets.greater_diameter_or_smaller_index(a,b, dim+1, true);
        };
        std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),cmp);
    }


    void compute_dim_0_pairs(std::vector<diameter_index_t>& all_edges,
                             std::vector<diameter_entry_t>& columns_to_reduce) {

        union_find dset(n);

        std::vector<diameter_index_t> edges;
        get_edges(edges,all_edges);
#ifdef DEBUG
        std::cerr << clear_line << std::flush;
        std::cout << "assembled " <<edges.size() << "/" << binomial_coeff(n, 2) <<  " edges" <<std::endl;
        std::cout << std::endl;
#endif
#ifdef PRINT_PERSISTENCE_PAIRS
        std::cerr << clear_line << std::flush;
        std::cout << "persistence intervals in dim 0:" << std::endl;
#endif

        static simplex_coboundary_enumerator cofacets(*this);

        auto cmp = [&]
                (const diameter_entry_t& a, const diameter_entry_t& b) {
            return cofacets.greater_diameter_or_smaller_index(a,b, 1, true);
        };
        std::sort(edges.rbegin(), edges.rend(),cmp);

        std::vector<index_t> vertices_of_edge(2);
        for (auto e : edges) {

            get_simplex_vertices(get_index(e), 1, n, vertices_of_edge.rbegin(), binomial_coeff);
            index_t u = dset.find(vertices_of_edge[0]), v = dset.find(vertices_of_edge[1]);

            if (u != v) {
#ifdef PRINT_PERSISTENCE_PAIRS
                print_pairs(true,0,get_diameter(e));
#endif
                dset.link(u, v);
            }
            else if (dim_max > 0) {
                columns_to_reduce.push_back(e);
            }
        }
        if (dim_max > 0) std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

#ifdef PRINT_PERSISTENCE_PAIRS
        for (index_t i = 0; i < n; ++i)
            if (dset.find(i) == i) std::cout << " [0, )" << std::endl;
#endif
    }


    template <typename Column> diameter_entry_t pop_pivot(Column& column) {
        diameter_entry_t pivot(-1);
#ifdef USE_COEFFICIENTS
        while (!column.empty()) {
			if (get_coefficient(pivot) == 0)
				pivot = column.top();
			else if (get_index(column.top()) != get_index(pivot))
				return pivot;
			else
				set_coefficient(pivot,(get_coefficient(pivot) + get_coefficient(column.top())) % modulus);
			column.pop();
		}
		return (get_coefficient(pivot) == 0) ? -1 : pivot;
#else
        while (!column.empty()) {
            pivot = column.top();
            column.pop();
            if (column.empty() || get_index(column.top()) != get_index(pivot)) return pivot;
            column.pop();
        }
        return -1;
#endif
    }


    template <typename Column> diameter_entry_t get_pivot(Column& column) {
        diameter_entry_t result = pop_pivot(column);
        if (get_index(result) != -1) column.push(result);
        return result;
    }


    template <typename Column>
    diameter_entry_t init_coboundary_and_get_pivot(const diameter_entry_t simplex,
                                                   Column& working_coboundary, const index_t& dim,
                                                   entry_hash_map& pivot_column_index) {

        static simplex_coboundary_enumerator cofacets(*this);
        cofacet_entries.clear();
        cofacets.set_simplex(simplex, dim, false, true);

        while (cofacets.has_next()) {
            diameter_entry_t cofacet = cofacets.next(false, true);
            if (!use_threshold ||is_lower(get_diameter(cofacet),threshold)) {
                cofacet_entries.push_back(cofacet);
            }
        }
        for (auto cofacet : cofacet_entries) working_coboundary.push(cofacet);
        return get_pivot(working_coboundary);
    }


    template <typename Column>
    void add_simplex_coboundary(const diameter_entry_t simplex, const index_t& dim,
                                std::stack<diameter_entry_t>& working_reduction_column, Column& working_coboundary) {

        static simplex_coboundary_enumerator cofacets(*this);
        working_reduction_column.push(simplex);
        cofacets.set_simplex(simplex, dim);
        while (cofacets.has_next()) {
            diameter_entry_t cofacet = cofacets.next();

            if (!use_threshold ||is_lower(get_diameter(cofacet),threshold)) working_coboundary.push(cofacet);
        }
    }


    template <typename Column>
    void add_coboundary(compressed_sparse_matrix<diameter_entry_t>& reduction_matrix,
                        const std::vector<diameter_entry_t>& columns_to_reduce,
                        const size_t index_column_to_add, const coefficient_t factor,
                        const size_t& dim, std::stack<diameter_entry_t>& working_reduction_column,
                        Column& working_coboundary) {
        diameter_entry_t column_to_add(columns_to_reduce[index_column_to_add], factor);
        add_simplex_coboundary(column_to_add, dim, working_reduction_column, working_coboundary);

        for (diameter_entry_t simplex : reduction_matrix.subrange(index_column_to_add)) {
            set_coefficient(simplex, get_coefficient(simplex) * factor % modulus);
            add_simplex_coboundary(simplex, dim, working_reduction_column, working_coboundary);
        }
    }


    void compute_pairs(const std::vector<diameter_entry_t>& columns_to_reduce,
                       entry_hash_map& pivot_column_index, const index_t dim) {

        bool threshold_reached= false;

        static simplex_coboundary_enumerator cofacets(*this);

#ifdef PRINT_PERSISTENCE_PAIRS
        std::cout << "persistence intervals in dim " << dim << ":" << std::endl;
#endif
        compressed_sparse_matrix<diameter_entry_t> reduction_matrix;

#ifdef INDICATE_PROGRESS
        std::chrono::steady_clock::time_point next = std::chrono::steady_clock::now() + time_step;
#endif
        for (size_t index_column_to_reduce = 0; index_column_to_reduce < columns_to_reduce.size();
             ++index_column_to_reduce) {

            threshold_reached= false;

            diameter_entry_t column_to_reduce(columns_to_reduce[index_column_to_reduce], 1);
            PT diameter = get_diameter(column_to_reduce);

            reduction_matrix.append_column();

            auto cmp = [&](const diameter_entry_t& a, const diameter_entry_t& b) {
                return cofacets.greater_diameter_or_smaller_index(a,b, dim+1, false);
            };

            std::priority_queue<diameter_entry_t, std::vector<diameter_entry_t>,
                    decltype(cmp)> working_coboundary(cmp);
            std::stack<diameter_entry_t> working_reduction_column;

            diameter_entry_t f, pivot = init_coboundary_and_get_pivot(
                    column_to_reduce, working_coboundary, dim, pivot_column_index);

            while (true) {
#ifdef INDICATE_PROGRESS
                if (std::chrono::steady_clock::now() > next) {
                    std::cerr << clear_line << "reducing column " << index_column_to_reduce + 1
                              << "/" << columns_to_reduce.size() << " (diameter [" << diameter
                              << "," << std::sqrt(CGAL::to_double(get_diameter(pivot)))/2
                              << "), index " << get_index(pivot) << " , steps " << step_count
                              << ")" << std::flush;
                    next = std::chrono::steady_clock::now() + time_step;
                }
#endif
                if (get_index(pivot) != -1) {
                    auto pair = pivot_column_index.find(get_entry(pivot));

                    if (pair != pivot_column_index.end()) {
                        entry_t other_pivot = pair->first;
                        index_t index_column_to_add = pair->second;
                        step_count++;
                        coefficient_t factor = modulus - get_coefficient(pivot) *
                                                         multiplicative_inverse[get_coefficient(other_pivot)] % modulus;

                        add_coboundary(reduction_matrix, columns_to_reduce, index_column_to_add,
                                       factor, dim, working_reduction_column, working_coboundary);
                        pivot = get_pivot(working_coboundary);

                    }
                    else if (get_index(f = get_zero_apparent_facet(pivot, dim + 1)) != -1) {
                        set_coefficient(f, modulus - get_coefficient(f));
                        step_count++;
                        add_simplex_coboundary(f, dim, working_reduction_column, working_coboundary);
                        pivot = get_pivot(working_coboundary);

                    }
                    else break;
                }
                else {
                    threshold_reached= true;
                    break;
                }
            }
            if (threshold_reached== true){
#ifdef PRINT_PERSISTENCE_PAIRS
#ifdef INDICATE_PROGRESS
                std::cerr << clear_line << std::flush;
#endif
                print_pairs(false,diameter,0);
#endif
            }
            else{
#ifdef PRINT_PERSISTENCE_PAIRS
                PT death = get_diameter(pivot);
#ifdef INDICATE_PROGRESS
                std::cerr << clear_line << std::flush;
#endif
                print_pairs(true,diameter,death);
#endif
                pivot_column_index.insert({ get_entry(pivot), index_column_to_reduce });

                while (true) {
                    if (working_reduction_column.empty()) break;
                    diameter_entry_t f = working_reduction_column.top();
                    working_reduction_column.pop();
                    reduction_matrix.push_back(f);
                }
            }

        }
#ifdef INDICATE_PROGRESS
        std::cerr << clear_line << std::flush;
#endif
    }

    void get_edges(std::vector<diameter_index_t>& edges, std::vector<diameter_index_t>& all_edges);

    void compute_barcodes() {
        std::vector<diameter_index_t> simplices;
        std::vector<diameter_entry_t> columns_to_reduce;

        compute_dim_0_pairs(simplices, columns_to_reduce);

        for (index_t dim = 1; dim <= dim_max; ++dim) {
            entry_hash_map pivot_column_index;
            pivot_column_index.reserve(columns_to_reduce.size());

            compute_pairs(columns_to_reduce, pivot_column_index, dim);

            if (dim < dim_max) {
                assemble_columns_to_reduce(simplices, columns_to_reduce, pivot_column_index, dim+1);
            }
        }
    }
};


template <> class cecher<distance_matrix>::simplex_coboundary_enumerator {
    index_t idx_below, idx_above, j, k;
    std::vector<index_t> vertices;
    diameter_entry_t simplex;
    const coefficient_t modulus;
    const distance_matrix& dist;
    const table& binomial_coeff;
    const cecher& parent;
    const mtx_PT& points_PT;
    const mtx_ET& points_ET;
    int dim;

    min_ball<PT> min_ball_PT;
    min_ball<ET> min_ball_ET;
    min_2_ball<PT> min_2_ball_PT;
    bool simplex_degenerate;

public:

    simplex_coboundary_enumerator(const cecher& _parent) : modulus(_parent.modulus), dist(_parent.dist),
                                                           binomial_coeff(_parent.binomial_coeff), parent(_parent),
                                                           points_PT(_parent.points_PT), points_ET(_parent.points_ET),
                                                           min_ball_PT(min_ball<PT>(points_PT,binomial_coeff,dist)),
                                                           min_ball_ET(min_ball<ET>(points_ET,binomial_coeff,dist)),
                                                           min_2_ball_PT(min_2_ball<PT>(binomial_coeff,dist)){}

    bool is_subset(const std::vector<index_t>& supp_a,  const std::vector<index_t>& supp_b) {
        int l=0;
        for (int k=0; k< supp_a.size(); k++) {
            for (int j=0; j< supp_b.size(); j++) if (supp_b[j]==supp_a[k]) l++;
        }
        return (l==supp_a.size());
    }

    void next_indices(int& i, int& j, const int& max){
        if(j<max) j++;
        else {
            i++;
            j=i;
        }
    }

    bool greater_diameter_or_smaller_index(const diameter_entry_t& a, const diameter_entry_t& b,
                                           const index_t dim, const bool is_assembly) {

        if (!is_assembly && get_support(a)==get_support(b)) return (get_index(a) < get_index(b));
        else {
#ifndef USE_RATIONALS
            try {
#endif
                if (get_diameter(a)!=get_diameter(b)) return (get_diameter(a)>get_diameter(b));
#ifndef USE_RATIONALS
                throw CGAL::Uncertain_conversion_exception("");
            }
            catch(CGAL::Uncertain_conversion_exception) {
                recomputation_count++;

#endif
                std::vector<index_t> supp_a(get_supp_dim(a)+1), supp_b(get_supp_dim(b)+1);
                get_simplex_vertices(get_supp_index(a), get_supp_dim(a), parent.n, supp_a.rbegin(), binomial_coeff);
                get_simplex_vertices(get_supp_index(b), get_supp_dim(b), parent.n, supp_b.rbegin(), binomial_coeff);

#ifndef USE_RATIONALS
                upper_tri_mtx<ET> CM_inv_a= min_ball_ET.circumsphere(supp_a);
                upper_tri_mtx<ET> CM_inv_b= min_ball_ET.circumsphere(supp_b);

                if (CM_inv_a[0][0]!=CM_inv_b[0][0]) return (CM_inv_a[0][0]<CM_inv_b[0][0]);
#endif
                if (dim==1 || (get_supp_dim(a)==1 && get_supp_dim(b)==1)) {
                    return (supp_a[0]<supp_b[0] || (supp_a[0]==supp_b[0] && supp_a[1]<supp_b[1]));
                }

                if (supp_a.size()<supp_b.size() && is_subset(supp_a,supp_b)) return false;
                if (supp_b.size()<supp_a.size() && is_subset(supp_b,supp_a)) return true;
#ifdef USE_RATIONALS
                upper_tri_mtx<ET> CM_inv_a= min_ball_ET.circumsphere(supp_a);
                upper_tri_mtx<ET> CM_inv_b= min_ball_ET.circumsphere(supp_b);
#endif
                int i1=0,i2=1,j1=0,j2=1;
                while (true) {

                    while (i1 < supp_a.size()+1 && CM_inv_a[i2][i1]==0) next_indices(i1,i2,supp_a.size());
                    while (j1 < supp_b.size()+1 && CM_inv_b[j2][j1]==0) next_indices(j1,j2,supp_b.size());

                    if (i1>0 && j1==0) return false;
                    if (j1==0 && j1>0) return true;
                    if (i1==supp_a.size()&& j1<supp_b.size()) return (CM_inv_b[j2][j1]>0);
                    if (i1<supp_a.size()&& j1==supp_b.size()) return (CM_inv_a[i2][i1]<0);
                    if (i1==supp_a.size()&& j1==supp_b.size()) return(get_index(a) < get_index(b));

                    if (i1==0 && supp_a[i2-1]!=supp_b[j2-1]) return (supp_a[i2-1]<supp_b[j2-1]);
                    if (i1>0 && (supp_a[i1-1]>supp_b[j1-1] || (supp_a[i1-1]==supp_b[j1-1] && supp_a[i2-1]>supp_b[j2-1]))) {
                        return (CM_inv_b[j2][j1]>0);
                    }
                    if (i1>0 && (supp_a[i1-1]<supp_b[j1-1] || (supp_a[i1-1]==supp_b[j1-1] && supp_a[i2-1]<supp_b[j2-1]))) {
                        return (CM_inv_a[i2][i1]<0);
                    }

                    if (i1==0 && CM_inv_a[i2][i1]!=CM_inv_b[j2][j1]) return (CM_inv_a[i2][i1]>CM_inv_b[j2][j1]);
                    if (i1>0 && CM_inv_a[i2][i1]!=CM_inv_b[j2][j1]) return (CM_inv_a[i2][i1]<CM_inv_b[j2][j1]);

                    next_indices(i1,i2,supp_a.size());
                    next_indices(j1,j2,supp_b.size());
                }
#ifndef USE_RATIONALS
            }
#endif
        }
    }


    void set_simplex(const diameter_entry_t _simplex, const index_t _dim,
                     const bool is_assembly=false, const bool init_coboundary=false) {
        dim =_dim;
        idx_below = get_index(_simplex);
        idx_above = 0;
        j = parent.n - 1;
        k = _dim + 1;
        simplex = _simplex;
        vertices.resize(_dim + 1);
        parent.get_simplex_vertices(get_index(_simplex), _dim, parent.n, vertices.rbegin(), binomial_coeff);
        simplex_degenerate=false;
        min_ball_ET.set(is_assembly,init_coboundary);
        min_ball_PT.set(is_assembly,init_coboundary);
        min_2_ball_PT.set(_simplex,init_coboundary);
#ifndef USE_RATIONALS
        if (dim>1 || (dim>0) && is_assembly && has_next(false)) {
            try {min_ball_PT.construct_simplex_cache(vertices);}
            catch(CGAL::Uncertain_conversion_exception) {simplex_degenerate=true;}
        }
#endif
    }

    bool has_next(bool all_cofacets = true) {
        return (j >= k && (all_cofacets || binomial_coeff(j, k) > idx_below));
    }

    diameter_entry_t next(const bool is_assembly=false ,const bool init_coboundary=false) {

        while ((binomial_coeff(j, k) <= idx_below)) {
            idx_below -= binomial_coeff(j, k);
            idx_above += binomial_coeff(j, k + 1);
            --j;
            --k;
            assert(k != -1);
        }

        index_t cofacet_index = idx_above + binomial_coeff(j, k + 1) + idx_below;
        coefficient_t cofacet_coeff =
                (k & 1 ? modulus - 1 : 1) * get_coefficient(simplex) % modulus;

        if (simplex_degenerate) return diameter_entry_t(min_ball_ET.construct(vertices,j--), cofacet_index, cofacet_coeff);
        try {
            if (dim>1 || (dim>0) && is_assembly) {
#ifdef USE_RATIONALS
                return diameter_entry_t(min_ball_PT.construct(vertices, j--), cofacet_index, cofacet_coeff);
#else
                return diameter_entry_t(min_ball_PT.construct_cofacet(vertices, j--), cofacet_index,  cofacet_coeff);
#endif
            }
            else return diameter_entry_t(min_2_ball_PT.construct(vertices, j--,cofacet_index),cofacet_index, cofacet_coeff); 
        }
        catch(CGAL::Uncertain_conversion_exception) {
            recomputation_count++;
            return diameter_entry_t(min_ball_ET.construct(vertices, j+1), cofacet_index, cofacet_coeff);
        }
    }
};


template <> void cecher<distance_matrix>::get_edges(std::vector<diameter_index_t>& edges,
                                                    std::vector<diameter_index_t>& all_edges) {

    std::vector<index_t> vertices(2);

#ifdef INDICATE_PROGRESS
    std::chrono::steady_clock::time_point next = std::chrono::steady_clock::now() + time_step;
#endif
    vertices[0]=dist.size()-2;
    vertices[1]=dist.size()-1;

    for (index_t index = binomial_coeff(n, 2); index-- > 0;) {
        PT length = dist(vertices[0], vertices[1]);

#ifdef INDICATE_PROGRESS
        if (std::chrono::steady_clock::now() > next) {
            std::cerr << clear_line << "assembling"
                      << " edges (processing " << binomial_coeff(n, 2)-index
                      << "/" << binomial_coeff(n, 2) << " simplices)" << std::flush;
            next = std::chrono::steady_clock::now() + time_step;
        }
#endif

        if (!use_threshold ||is_lower(length,threshold)) {

            all_edges.push_back({ length, index });
            bool is_apparent=false;
            for (index_t l = n-1; l >=0; --l) {
                if (l!=vertices[0] && l!=vertices[1]){
                    try{
                        PT orientation= dist(vertices[0],l) + dist(vertices[1],l) - length;
                        if (orientation < 0 || orientation ==0 && l>vertices[0]) {
                            is_apparent=true;
                            break;
                        }
                    }
                    catch(CGAL::Uncertain_conversion_exception) {
                        recomputation_count++;
                        ET orientation= -compute_dist(points_ET[vertices[0]],points_ET[vertices[1]]);
                        for (int i=0; i<2; i++) orientation+=compute_dist(points_ET[vertices[i]],points_ET[l]);
                        if (orientation < 0 || orientation ==0 && l>vertices[0]) {
                            is_apparent=true;
                            break;
                        }
                    }
                }
            }
            if (is_apparent==false) edges.push_back({ length, index });
        }
        if (vertices[0]==0) {
            --vertices[1];
            vertices[0]=vertices[1]-1;
        }
        else --vertices[0];
    }
}

template <typename T> class read_convert {
public:
    static T convert_from_rational(const ET cofacet_diameter);

    static std::vector<std::vector<T>> read_points(std::istream& input_stream) {
        std::vector<std::vector<T>> points;
        std::string line;
        ET value;
        while (std::getline(input_stream, line)) {
            std::vector<T> point;
            std::istringstream s(line);
            while (s >> value) {
                point.push_back(convert_from_rational(value));
                s.ignore();
            }
            if (!point.empty()) points.push_back(point);
            assert(point.size() == points.front().size());
        }
        return points;
    }
};

template<>ET read_convert<ET>::convert_from_rational(const ET val){
    return val;
}
template<>CGAL::Interval_nt<> read_convert<CGAL::Interval_nt<>>::convert_from_rational(const ET val){
    return CGAL::to_interval(val);
}
template<>double read_convert<double>::convert_from_rational(const ET val){
    return CGAL::to_double(val);
}
template<>float read_convert<float>::convert_from_rational(const ET val){
    return static_cast<float>(CGAL::to_double(val));
}


euclidean_distance_matrix read_point_cloud(std::istream& input_stream) {

    std::vector<std::vector<PT>> points=read_convert<PT>::read_points(input_stream);
    euclidean_distance_matrix eucl_dist(std::move(points));

    index_t n = eucl_dist.size();
    point_dim =eucl_dist.points_PT.front().size();
    std::cout << "point cloud with " << n << " points in dimension "
              <<  point_dim << std::endl;

    return eucl_dist;
}

distance_matrix read_file(std::istream& input_stream) {
    return read_point_cloud(input_stream);
}

void print_usage_and_exit(int exit_code) {
    std::cerr
            << "Usage: "
            << "Cecher "
            << "[options] [filename]" << std::endl
            << std::endl
            << "Options:" << std::endl
            << std::endl
            << "  --help           print this screen" << std::endl
            << "  --dim <k>        compute Cech persistent homology up to dimension k" << std::endl
            << "  --threshold <t>  compute Cech complexes up to radius t (decimal)" << std::endl
            #ifdef USE_COEFFICIENTS
            << "  --modulus <p>    compute homology with coefficients in the prime field Z/pZ"<< std::endl
            #endif
            << std::endl;
    exit(exit_code);
}

int main(int argc, char** argv) {
#ifdef DEBUG
    auto start = std::chrono::steady_clock::now();
#endif
    index_t dim_max = 1;
    PT threshold = std::numeric_limits<PT>::max();
    PT ratio = 1;
    coefficient_t modulus = 2;

    for (index_t i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--help") {
            print_usage_and_exit(0);
        }
        else if (arg == "--dim") {

            std::string parameter = std::string(argv[++i]);
            size_t next_pos;
            dim_max = std::stol(parameter, &next_pos);

            if (next_pos != parameter.size()) print_usage_and_exit(-1);
        }
        else if (arg == "--threshold") {
            use_threshold = true;
            std::string parameter = std::string(argv[++i]);
            size_t next_pos;

            threshold = std::stof(parameter, &next_pos);
            threshold = threshold * threshold *4;

            if (next_pos != parameter.size()) print_usage_and_exit(-1);

#ifdef USE_COEFFICIENTS
            }
		else if (arg == "--modulus") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			modulus = std::stol(parameter, &next_pos);
			if (next_pos != parameter.size() || !is_prime(modulus)) print_usage_and_exit(-1);
#endif
        }
        else {
            if (filename) { print_usage_and_exit(-1); }
            filename = argv[i];
        }
    }

    std::ifstream file_stream(filename);
    if (filename && file_stream.fail()) {
        std::cerr << "couldn't open file " << filename << std::endl;
        exit(-1);
    }
    distance_matrix dist = read_file(filename ? file_stream : std::cin);

    std::ifstream file_stream_ET(filename);
    mtx_ET points_ET=read_convert<ET>::read_points(filename ? file_stream_ET : std::cin);
    std::ifstream file_stream_PT(filename);
    mtx_PT points_PT=read_convert<PT>::read_points(filename ? file_stream_PT : std::cin);

    cecher<distance_matrix>(std::move(dist), dim_max,
                            threshold,ratio, modulus, std::move(points_PT), std::move(points_ET)).compute_barcodes();

#ifdef DEBUG
    auto end = std::chrono::steady_clock::now();
    double elapsed_seconds = std::chrono::duration_cast< std::chrono::duration<double> >(end - start).count();
    std::cout << std::endl;
    std::cout << "### process finished in: " << elapsed_seconds << "s ###" <<  '\n';
    std::cout << "recomputations: " << recomputation_count <<  '\n';
    std::cout << "steps: " << step_count <<  '\n';
#endif
    exit(0);
}
