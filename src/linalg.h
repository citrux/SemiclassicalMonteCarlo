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

inline vec2 operator-(vec2 v) { return {-v.x, -v.y}; }

inline vec2 operator+(vec2 a, vec2 b) { return {a.x + b.x, a.y + b.y}; }

inline vec2 operator-(vec2 a, vec2 b) { return a + (-b); }

inline vec2 operator*(vec2 v, double scale) {
    return {scale * v.x, scale * v.y};
}

inline vec2 operator*(double scale, vec2 v) { return v * scale; }

inline vec2 operator/(vec2 v, double scale) { return v * (1. / scale); }

inline double dot(vec2 a, vec2 b) { return a.x * b.x + a.y * b.y; }

inline double cross(vec2 a, vec2 b) { return a.x * b.y - a.y * b.x; }

inline double len(vec2 v) { return sqrt(dot(v, v)); }

inline vec2 ort(vec2 v) { return v / len(v); }

inline vec2 operator-(Point end, Point start) {
    return {end.x - start.x, end.y - start.y};
}

inline Point operator+(Point start, vec2 shift) {
    return {start.x + shift.x, start.y + shift.y};
}

inline Point operator-(Point start, vec2 shift) { return start + (-shift); }