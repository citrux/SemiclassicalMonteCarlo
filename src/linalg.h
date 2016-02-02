/*
    *
    * Структура, хранящая компоненты импульса
    *
*/
struct Point {
    double x;
    double y;
};

struct vec2 {
    double x;
    double y;
};

inline vec2 operator-(const vec2 & v) { return {-v.x, -v.y}; }

inline vec2 operator+(const vec2 & a, const vec2 & b) {
    return {a.x + b.x, a.y + b.y};
}

inline vec2 operator-(const vec2 & a, const vec2 & b) { return a + (-b); }

inline vec2 operator*(const vec2 v, double scale) {
    return {scale * v.x, scale * v.y};
}

inline vec2 operator*(const double scale, const vec2 & v) { return v * scale; }

inline vec2 operator/(const vec2 & v, double scale) { return v * (1. / scale); }

inline double dot(const vec2 & a, const vec2 & b) {
    return a.x * b.x + a.y * b.y;
}

inline double cross(const vec2 & a, const vec2 & b) {
    return a.x * b.y - a.y * b.x;
}

inline double len(const vec2 & v) { return sqrt(dot(v, v)); }

inline vec2 ort(const vec2 & v) { return v / len(v); }

inline vec2 operator-(const Point & end, const Point & start) {
    return {end.x - start.x, end.y - start.y};
}

inline Point & operator+=(Point & start, const vec2 & shift) {
    start.x += shift.x;
    start.y += shift.y;
    return start;
}
inline Point operator+(Point start, const vec2 & shift) {
    return start += shift;
}
inline Point operator-(Point start, const vec2 & shift) {
    return start += (-shift);
}