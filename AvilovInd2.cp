#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
//Avilov Alex - Backward Ray Tracing
using namespace::std;

#define PI 3.14
#define MAX_RAY_DEPTH 5

//Класс вектора в трехмерном пространстве
class Vec3
{
public:
    float x, y, z;
    Vec3() : x(0.0), y(0.0), z(0.0) {}
    Vec3(float xx) : x(xx), y(xx), z(xx) {}
    Vec3(float xx, float yy, float zz) : x(xx), y(yy), z(zz) {}
    Vec3& normalize()
    {
        float nor2 = length2();
        if (nor2 > 0) {
            float invNor = 1 / sqrt(nor2);
            x *= invNor, y *= invNor, z *= invNor;
        }
        return *this;
    }
    Vec3 operator * (const float &f) const { return Vec3(x * f, y * f, z * f); }
    Vec3 operator * (const Vec3 &v) const { return Vec3(x * v.x, y * v.y, z * v.z); }
    float dot(const Vec3 &v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3 operator - (const Vec3 &v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator + (const Vec3 &v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3& operator += (const Vec3 &v) { x += v.x, y += v.y, z += v.z; return *this; }
    Vec3& operator *= (const Vec3 &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
    Vec3 operator - () const { return Vec3(-x, -y, -z); }
    float length2() const { return x * x + y * y + z * z; }
    float length() const { return sqrt(length2()); }
    friend ostream & operator << (ostream &os, const Vec3 &v)
    {
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
};

class Sphere
{
public:
    Vec3 center;
    float radius, radius2;
    //Цвет поверхности
    Vec3 surfaceColor;
    //Цвет излучения
    Vec3 emissionColor;
    //Прозрачность
    float transparency;
    //Отражение
    float reflection;
    Sphere(
        const Vec3 &c,
        const float &r,
        const Vec3 &sc,
        const float &refl = 0,
        const float &transp = 0.5,
        const Vec3 &ec = 0) :
        center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
        transparency(transp), reflection(refl) {}
//Пересечение
bool intersect(const Vec3 &rayorig, const Vec3 &raydir, float &t0, float &t1) const
    {
        Vec3 v = center - rayorig;
        float b = v.dot(raydir); //b
        if (b < 0) return false;  //<0 -> пересечения нет
        float c = v.dot(v) - b * b;//нужно чтобы проверить лежит ли начало луча в сфере
        if (c > radius2) return false;
        float thc = sqrt(radius2 - c); //
        t0 = b - thc; //точка пересечения
        t1 = b + thc; //еще точка пересечения

        return true;
    }
};

float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}
//Берем луч, определенный его направлением и исходной точкой -> находим пересечения с объектами сцены
//-> если пересекает, то находим точки пересечения, нормаль в этих точках и считаем там цвет
Vec3 trace(const Vec3 &rayorig, const Vec3 &raydir, const vector<Sphere> &objects, const int &depth) {
    float tnear = 1000000;
    const Sphere* sphere = NULL;
    //сначала найти пересечение луча со сферой
    for (unsigned i = 0; i < objects.size(); ++i) {
        float t0 = 1, t1 = 1;
        if (objects[i].intersect(rayorig, raydir, t0, t1)) {
            if (t0 < 0) t0 = t1;
            if (t0 < tnear) {
                tnear = t0;
                sphere = &objects[i];
            }
        }
    }
    //если пересечения нет, то цвет фона
    if (!sphere) return Vec3(2);
    Vec3 surfaceColor = 0; //цвет пересеченного объекта
    Vec3 phit = rayorig + raydir * tnear; //точка пересечения
    Vec3 nhit = phit - sphere->center; //нормаль
    nhit.normalize();
    bool inside = false;
    // если нормаль и направление луча не противоположны
    if (raydir.dot(nhit) > 0) //поменяем направление нормали в точке пересечения, при этом мы находимся внутри сферы, поэтому inside = true
        nhit = -nhit;inside = true;
    //если глубина позволяет, то просчитываем прозрачность и отражение
    if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
        float fresneleffect = mix(pow(1 - (-raydir.dot(nhit)), 3), 1, 0.1); //какой-то эффект отражения
        // считаем отражение
        Vec3 refldir = raydir - nhit * 2 * raydir.dot(nhit); //направление отраженного луча (первичный луч - 2 * вектор нормали в точке пересечения * скалярное произведение)
        refldir.normalize();
        Vec3 reflection = trace(phit + nhit, refldir, objects, depth + 1); //строим отраженный вектор
        Vec3 refraction = 0;
        //прозрачность
        if (sphere->transparency) {
            float ior = 1.1; // коэффициент рефракции (второй просто равен 1)
            float eta = (inside) ? ior : 1 / ior; //частное коэффициентов рефракции для формулы. Зависит от среды(снаружи или внутри сферы)
            float k = -nhit.dot(raydir);//скалярное произведение для формулы
            float cosi = 1 - eta * eta * (1 - k * k);//угол между нормалью и преломленным лучом
            Vec3 refrdir = raydir * eta + nhit * (eta *  k - sqrt(cosi)); //направление по формуле направления
            refrdir.normalize();
            refraction = trace(phit - nhit, refrdir, objects, depth + 1); //строим вектор в направлении преломления прозрачной поверхности
        }
        // в результате смесь отражения и прозрачности
        surfaceColor = (
            reflection * fresneleffect + refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surfaceColor; //интенсивность отражения, прозрачности и исходный цвет объекта
    }
    else {
        for (unsigned i = 0; i < objects.size(); ++i) {
            if (objects[i].emissionColor.x > 0) {
                //свет
                Vec3 transmission = 1;
                Vec3 lightDirection = objects[i].center - phit;
                lightDirection.normalize();
                for (unsigned j = 0; j < objects.size(); ++j) {
                    if (i != j) {
                        float t0, t1;
                        if (objects[j].intersect(phit + nhit, lightDirection, t0, t1)) {
                            transmission = 0;
                            break;
                        }
                    }
                }
                surfaceColor += sphere->surfaceColor * transmission * max(float(0), nhit.dot(lightDirection)) * objects[i].emissionColor;
            }
        }
    }

    return surfaceColor + sphere->emissionColor;
}
//для каждого пиксела считаем цвет(сферы или фона)
 void render(const vector<Sphere> &objects) {
    unsigned width = 640, height = 480;
    Vec3 *image = new Vec3[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30; // расстояние камеры от сетки экрана
    float angle = tan(PI * 0.5 * fov / 180);
    for (unsigned y = 0; y < height; ++y) {
        for (unsigned x = 0; x < width; ++x, ++pixel) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * (width / float(height));
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3 raydir(xx, yy, -1);
            raydir.normalize();
            *pixel = trace(Vec3(0), raydir, objects, 0);
        }
    }
    ofstream ofs("./macButton.ppm", ios::out | ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (unsigned i = 0; i < width * height; ++i) {
        ofs << (unsigned char)(min(float(1), image[i].x) * 255) <<
               (unsigned char)(min(float(1), image[i].y) * 255) <<
               (unsigned char)(min(float(1), image[i].z) * 255);
    }
    ofs.close();
    delete [] image;
}
 int main(int argc, char **argv)
{
    vector<Sphere> objects;
    objects.push_back(Sphere(Vec3(-7.0, 2, -50), 6, Vec3(1.00, 0.384, 0.35)));
    objects.push_back(Sphere(Vec3(0.0, 0, -20), 2, Vec3(1, 0.74, 0.18)));
    objects.push_back(Sphere(Vec3(3.0, -1, -30), 3, Vec3(0.156, 0.807, 0.25)));

    //cцена - это большая сфера снизу
    objects.push_back(Sphere(Vec3(0.0, -1000004, -20), 1000000, Vec3(0.20, 0.20, 0.20), 0, 0.0));
    //свет - это маленькаы сфера сверху (не отражающая и не прозрачная)
    objects.push_back(Sphere(Vec3(0.0, 30, -30), 1, Vec3(0.00, 0.00, 0.00), 0, 0.0, Vec3(3)));
    render(objects);

    return 0;
}
