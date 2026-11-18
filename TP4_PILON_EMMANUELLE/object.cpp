#include "object.hpp"

#include <cmath>
#include <cfloat>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>

// ---------------------------------------------------------
// 1. TRANSFORMATION GLOBALE (Object::intersect)
// ---------------------------------------------------------
// Cette fonction sert d'intermédiaire. Le rayon arrive en "coordonnées monde".
// Pour simplifier les maths, on transforme le rayon en "coordonnées objet" (local),
// on teste l'intersection, puis on retransforme le résultat en "coordonnées monde".
bool Object::intersect(Ray ray, Intersection &hit) const 
{
    // [EXPL] On s'assure que w=1 pour la position (point) et w=0 pour la direction (vecteur)
    ray.origin[3] = 1;
    ray.direction[3] = 0;

    // [EXPL] Transformation INVERSE : On amène le rayon dans l'espace local de l'objet.
    // Ex: Si l'objet est tourné de 90°, on tourne le rayon de -90°.
    Ray local_ray(i_transform * ray.origin, i_transform * ray.direction);
    
    // save best global depth
    double previous_best_global_depth = hit.depth;
    
    // reset for local test
    hit.depth = DBL_MAX;  

    if (localIntersect(local_ray, hit)) 
    {
        hit.position[3] = 1;
        hit.normal[3] = 0;
        
        // [EXPL] Transformation DIRECTE : On ramène le point d'intersection et la normale dans le monde réel.
        Vector global_position = transform * hit.position;
        // [EXPL] Attention : les normales se transforment avec la transposée de l'inverse (n_transform).
        Vector global_normal = (n_transform * hit.normal).normalized();
        
        Vector toIntersection = global_position - ray.origin;
        double global_depth = toIntersection.length();
        
        // [EXPL] Si cette intersection est plus proche que la précédente trouvée, on la garde.
        if (global_depth > 0 && global_depth < previous_best_global_depth) {
            
            hit.depth = global_depth;
            hit.position = global_position;
            hit.normal = global_normal;
            return true;
        }
    }

    // restore depth if fail
    hit.depth = previous_best_global_depth;
    return false;
}

bool Sphere::localIntersect(Ray const &ray, Intersection &hit) const 
{
    Vector direction = ray.direction.normalized();
    
    Vector oc = ray.origin;
    
    // quadratic coeffs
    double a = direction.dot(direction); 
    double b = 2.0 * oc.dot(direction);
    double c = oc.dot(oc) - radius * radius;
    
    double discriminant = b * b - 4 * a * c;
    
    if (discriminant < 0) {
        return false;
    }
    
    double sqrtDiscriminant = sqrt(discriminant);
    double t1 = (-b - sqrtDiscriminant) / (2.0 * a);
    double t2 = (-b + sqrtDiscriminant) / (2.0 * a);
    
    // pick closest valid t
    double t = t1;
    if (t1 < 0) {
        t = t2;
    }
    
    if (t < 0 || t >= hit.depth) {
        return false;
    }
    
    hit.depth = t;
    hit.position = ray.origin + t * direction;
    hit.normal = hit.position.normalized();
    
    return true;
}


bool Plane::localIntersect(Ray const &ray, Intersection &hit) const
{
    Vector direction = ray.direction.normalized();
    Vector normal(0, 0, 1);
    
    double denom = normal.dot(direction);
    
    // check parallel
    if (ABS_FLOAT(denom) < 1e-6) {
        return false; 
    }
    
    double t = -(ray.origin.dot(normal)) / denom;
    
    if (t < 0 || t >= hit.depth) {
        return false;
    }
    
    hit.depth = t;
    hit.position = ray.origin + t * direction;
    
    if (denom > 0) {
        hit.normal = -normal;
    } else {
        hit.normal = normal;
    }
    
    return true;
}

bool Conic::localIntersect(Ray const &ray, Intersection &hit) const {
    Vector direction = ray.direction.normalized();
    Vector origin = ray.origin;
    
    // [EXPL] Calcul de la pente (m) et de l'ordonnée à l'origine (b) du rayon du cône
    // r(z) = m * z + b
    double m = (radius2 - radius1) / (zMax - zMin);
    double b = radius1 - m * zMin;
    
    // Raccourcis pour lisibilité
    double ox = origin[0]; double oy = origin[1]; double oz = origin[2];
    double dx = direction[0]; double dy = direction[1]; double dz = direction[2];
    
    // [EXPL] Coefficients de l'équation quadratique At^2 + Bt + C = 0
    // Dérivés de x^2 + y^2 = (mz + b)^2 en substituant x,y,z par le rayon.
    double a = dx*dx + dy*dy - m*m*dz*dz;
    double b_coef = 2*(ox*dx + oy*dy) - 2*m*(m*oz + b)*dz;
    double c = ox*ox + oy*oy - (m*oz + b)*(m*oz + b);
    
    double discriminant = b_coef*b_coef - 4*a*c;
    
    if (discriminant < 0) return false; // Pas d'intersection (le rayon passe à côté)
    
    double sqrtDisc = sqrt(discriminant);
    double t1 = (-b_coef - sqrtDisc) / (2*a);
    double t2 = (-b_coef + sqrtDisc) / (2*a);
    
    double t = -1;
    
    // find closest valid intersection in range 
    // [EXPL] On teste les deux solutions (entrée et sortie du cône)
    for (int i = 0; i < 2; i++) {
        double t_test = (i == 0) ? t1 : t2;
        if (t_test < 0 || t_test >= hit.depth) continue; // Derrière la caméra ou trop loin
        
        // [EXPL] Vérification de la hauteur (Clipping)
        // Le cône mathématique est infini, on le coupe entre zMin et zMax.
        Vector p = origin + t_test * direction;
        double z = p[2];
        if (z < zMin || z > zMax) continue;
        
        // On garde le t le plus proche
        if (t < 0 || t_test < t) t = t_test;
    }
    
    if (t < 0) return false; // Aucune intersection valide sur la surface latérale
    
    // [EXPL] Si on arrive ici, on a touché le flanc du cône.
    hit.depth = t;
    hit.position = origin + t * direction;
    
    // [EXPL] Calcul de la normale pour un cône
    // Normale = vecteur (x, y, -r*m). C'est perpendiculaire à la pente.
    double z = hit.position[2];
    double r = m * z + b;
    Vector normal(hit.position[0], hit.position[1], -r * m);
    normal.normalize();
    
    // [EXPL] On s'assure que la normale pointe vers l'extérieur
    if (normal.dot(direction) > 0) normal = -normal;
    hit.normal = normal;
    
    // [EXPL] GESTION DES "CAPS" (Couvercles)
    // Ici, tu vérifies si le rayon tape les disques plats aux extrémités.
    if (ABS_FLOAT(dz) > 1e-6) { // Si le rayon n'est pas parallèle aux couvercles
        double t_cap1 = (zMin - oz) / dz; // Test pour zMin (devrait aussi tester zMax théoriquement)

        // Si le couvercle est plus proche que le flanc (t) ou si on n'a pas touché le flanc (t<0)
        if (t_cap1 > 0 && t_cap1 < hit.depth) {
            Vector p_cap = origin + t_cap1 * direction;

            // Est-ce qu'on est DANS le cercle du couvercle ? (x^2 + y^2 <= r^2)
            if (p_cap[0]*p_cap[0] + p_cap[1]*p_cap[1] <= radius1*radius1) {
                if (t_cap1 < t || t < 0) { 
                     hit.depth = t_cap1; hit.position = p_cap;
                     // La normale d'un couvercle est simplement l'axe Z (0,0,1) ou (0,0,-1)
                     hit.normal = Vector(0,0, (dz>0)?-1:1);
                     if (dz > 0) hit.normal = -Vector(0,0,1); else hit.normal = Vector(0,0,1);
                     return true;
                }
            }
        }
    }
    return true;
}

bool Mesh::localIntersect(Ray const &ray, Intersection &hit) const
{
    // aabb test
    double tNear = -DBL_MAX, tFar = DBL_MAX;
    for (int i = 0; i < 3; i++) {
        if (ABS_FLOAT(ray.direction[i]) < 1e-6) {
            if (ray.origin[i] < bboxMin[i] || ray.origin[i] > bboxMax[i]) return false;
        }
        else {
            double t1 = (bboxMin[i] - ray.origin[i]) / ray.direction[i];
            double t2 = (bboxMax[i] - ray.origin[i]) / ray.direction[i];
            if (t1 > t2) std::swap(t1, t2);
            if (t1 > tNear) tNear = t1;
            if (t2 < tFar) tFar = t2;
            if (tNear > tFar || tFar < 0) return false;
        }
    }

    bool isHit = false;
    for (size_t tri_i = 0; tri_i < triangles.size(); tri_i++) {
        if (intersectTriangle(ray, triangles[tri_i], hit)) {
            isHit = true;
        }
    }
    return isHit;
}

double Mesh::implicitLineEquation(double p_x, double p_y,
    double e1_x, double e1_y,
    double e2_x, double e2_y) const
{
    return (e2_y - e1_y)*(p_x - e1_x) - (e2_x - e1_x)*(p_y - e1_y);
}

bool Mesh::intersectTriangle(Ray const &ray,
    Triangle const &tri,
    Intersection &hit) const
{
    Vector const &p0 = positions[tri[0].pi];
    Vector const &p1 = positions[tri[1].pi];
    Vector const &p2 = positions[tri[2].pi];

    Vector direction = ray.direction.normalized();
    
    // moller trumbore algo
    Vector edge1 = p1 - p0;
    Vector edge2 = p2 - p0;
    
    Vector h = direction.cross(edge2);
    double a = edge1.dot(h);
    
    if (ABS_FLOAT(a) < 1e-6) return false;
    
    double f = 1.0 / a;
    Vector s = ray.origin - p0;
    
    // barycentric u
    double u = f * s.dot(h);
    
    if (u < 0.0 || u > 1.0) return false;
    
    Vector q = s.cross(edge1);
    
    // barycentric v
    double v = f * direction.dot(q);
    
    if (v < 0.0 || u + v > 1.0) return false;
    
    double t = f * edge2.dot(q);
    
    if (t < 0 || t >= hit.depth) return false;
    
    hit.depth = t;
    hit.position = ray.origin + t * direction;
    
    Vector normal = edge1.cross(edge2).normalized();
    if (normal.dot(direction) > 0) normal = -normal;
    hit.normal = normal;
    
    return true;
}