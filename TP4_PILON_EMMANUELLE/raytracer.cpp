#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>

#include "raytracer.hpp"
#include "image.hpp"

// Boucle principale
void Raytracer::render(const char *filename, const char *depth_filename,
                       Scene const &scene)
{
    Image colorImage(scene.resolution[0], scene.resolution[1]);
    Image depthImage(scene.resolution[0], scene.resolution[1]);
    
    double *zBuffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++) {
        zBuffer[i] = DBL_MAX;
    }
    
    double fovyRad = deg2rad(scene.camera.fovy);
    double viewPlaneHeight = 2.0 * tan(fovyRad / 2.0);
    double viewPlaneWidth = viewPlaneHeight * scene.camera.aspect;

    // cam coords
    // [EXPL] Configuration de la "caméra virtuelle".
    // On calcule la taille du plan image virtuel placé devant la caméra.
    // u, v, w sont les vecteurs qui définissent l'orientation de la caméra (repère local).
    Vector w = (scene.camera.position - scene.camera.center).normalized(); 
    Vector u = scene.camera.up.cross(w).normalized(); 
    Vector v = w.cross(u);
    // Coin bas-gauche de l'écran virtuel
    Vector viewPlaneBottomLeft = scene.camera.position - w - (viewPlaneWidth / 2.0) * u - (viewPlaneHeight / 2.0) * v;

    // [EXPL] Double boucle pour parcourir chaque pixel de l'image (X, Y)
    for(int y = 0; y < scene.resolution[1]; y++) {
        for(int x = 0; x < scene.resolution[0]; x++) {

            // [EXPL] Génération du RAYON PRIMAIRE
            // On trace une ligne qui part de l'œil et passe par le centre du pixel actuel.
            Ray ray;
            // calcul de direction u rayon
            if (scene.objects.empty())
            {
                ray = Ray(scene.camera.position, (Vector(-320, -320, 640) + Vector(x + 0.5, y + 0.5, 0) - scene.camera.position).normalized());
            }
            else
            {
                // calc pixel pos on viewplane
                double pixelU = (x + 0.5) / scene.resolution[0];
                double pixelV = (y + 0.5) / scene.resolution[1];

                Vector pixelPosition = viewPlaneBottomLeft + pixelU * viewPlaneWidth * u + pixelV * viewPlaneHeight * v;

                Vector direction = (pixelPosition - scene.camera.position).normalized();

                ray = Ray(scene.camera.position, direction);
            }

            // [EXPL] Lancement du rayon dans la scène !
            int rayDepth = 0; // Compteur pour limiter les récursions (reflets)
            Vector color(0,0,0,0);

            double depth = scene.camera.zFar;

            trace(ray, rayDepth, scene, color, depth);

            // depth test and zbuffer update
            // [EXPL] Z-Buffer : Si ce point est le plus proche vu jusqu'à présent, on dessine.
            if(depth >= scene.camera.zNear && depth <= scene.camera.zFar && 
                depth < zBuffer[x + y*scene.resolution[0]]) {

                zBuffer[x + y*scene.resolution[0]] = depth;

                colorImage.setPixel(x, y, color);
                depthImage.setPixel(x, y, (depth-scene.camera.zNear) / 
                                        (scene.camera.zFar-scene.camera.zNear));
            }
        }

        if (y % 100 == 0)
        {
            printf("Row %d pixels finished.\n", y);
        }
    }

    colorImage.writePNG(filename);
    depthImage.writePNG(depth_filename);

    printf("Ray tracing finished with images saved.\n");

    delete[] zBuffer;
}

// Trouve l'objet le plus proche touché par le rayon.
bool Raytracer::trace(Ray const &ray, 
                 int &rayDepth,
                 Scene const &scene,
                 Vector &outColor, double &depth)
{
    rayDepth++; // On incrémente la profondeur de récursion
    bool foundIntersection = false;
    
    // scene vide par defaut
    if (scene.objects.empty())
    {
        // default scene cube
        double x = 1280 / ray.direction[2] * ray.direction[0] + ray.origin[0];
        double y = 1280 / ray.direction[2] * ray.direction[1] + ray.origin[1];
        if ((x <= 160) && (x >= -160) && (y <= 160) && (y >= -160))
        {
            Material m; 
            m.emission = Vector(16.0, 0, 0); 
            m.reflect = 0;
            Intersection intersection;
            outColor = shade(ray, rayDepth, intersection, m, scene);
            depth = 1280;
            foundIntersection = true;
        }
    }
    else
    {
        Intersection closestHit;
        closestHit.depth = depth;
        Object* closestObject = nullptr;
        Material closestMaterial;

        // [EXPL] On teste l'intersection avec TOUS les objets de la scène.
        // C'est l'étape la plus coûteuse du Ray Tracing naïf.
        // find closest object
        for (size_t i = 0; i < scene.objects.size(); i++) {
            Intersection hit;
            hit.depth = closestHit.depth; 
            
            if (scene.objects[i]->intersect(ray, hit)) {
                // On garde seulement l'objet le plus proche (le plus petit 'depth')
                if (hit.depth < closestHit.depth) {
                    closestHit = hit;
                    closestObject = scene.objects[i];
                    closestMaterial = scene.objects[i]->material;
                }
            }
        }

        // [EXPL] Si on a touché quelque chose, on calcule sa couleur (Shading)
        if (closestObject != nullptr) {
            depth = closestHit.depth;
            outColor = shade(ray, rayDepth, closestHit, closestMaterial, scene);
            foundIntersection = true;
        }
    }

    rayDepth--;

    return foundIntersection;
}

// Modèle de Phong (Ambiant + Diffus + Spéculaire) + Ombres + Reflets + Réfraction
Vector Raytracer::shade(Ray const &ray,
                 int &rayDepth,
                 Intersection const &intersection,
                 Material const &material,
                 Scene const &scene)
{
    Vector diffuse(0);
    Vector ambient(0);
    Vector specular(0);
    
    // [EXPL] Pour chaque lumière dans la scène...
    for (auto lightIter = scene.lights.begin(); lightIter != scene.lights.end(); lightIter++)
    {
        const PointLight &light = *lightIter;
        
        // Ambiante : Lumière de base présente partout
        ambient += light.ambient * material.ambient;
        
        Vector toLight = (light.position - intersection.position);
        double distanceToLight = toLight.length();
        Vector L = toLight.normalized();
        
        // shadow ray offset
        // [EXPL] OMBRES (Shadow Ray)
        // On lance un rayon depuis le point d'impact vers la lumière.
        double epsilon = 1e-3; // Petit décalage pour ne pas s'auto-intersecter (acné)
        Vector shadowOrigin = intersection.position + epsilon * intersection.normal;
        Ray shadowRay(shadowOrigin, L);
        
        bool inShadow = false;
        Intersection shadowHit;
        shadowHit.depth = distanceToLight;
        // Si ce rayon touche un objet ENTRE le point et la lumière...
        for (size_t i = 0; i < scene.objects.size(); i++) {
            if (scene.objects[i]->intersect(shadowRay, shadowHit)) { // ...alors on est à l'ombre.
                inShadow = true;
                break;
            }
        }
        
        if (inShadow) { // Si à l'ombre, on ignore diffuse/spéculaire pour cette lumière.
            continue;
        }
        
        // [EXPL] PHONG : Calcul de l'éclairage Diffus et Spéculaire
        Vector N = intersection.normal;
        Vector V = (ray.origin - intersection.position).normalized();
        
        double a0 = light.attenuation[0];
        double a1 = light.attenuation[1];
        double a2 = light.attenuation[2];
        double attenuation = 1.0 / (a0 + a1 * distanceToLight + a2 * distanceToLight * distanceToLight);
        
        // phong lighting
        double NdotL = N.dot(L);
        if (NdotL > 0) {
            // Diffus : dépend de l'angle entre la normale et la lumière (Loi de Lambert)
            diffuse += attenuation * light.diffuse * material.diffuse * NdotL;
            
            // Spéculaire : le petit point brillant (reflet de la source lumineuse)
            Vector H = (L + V).normalized(); // Vecteur Blinn-Phong
            double NdotH = N.dot(H);
            if (NdotH > 0) {
                double specularFactor = pow(NdotH, material.shininess);
                specular += attenuation * light.specular * material.specular * specularFactor;
            }
        }
    }
    
    // [EXPL] RÉFLEXION (Miroir)
    Vector reflectedLight(0);
    // mirror reflection
    if ((material.reflect > 1e-6) && (rayDepth < MAX_RAY_RECURSION))
    {
        // Calcul du vecteur réfléchi R
        Vector D = ray.direction.normalized();
        Vector N = intersection.normal;
        Vector R = D - 2.0 * D.dot(N) * N;
        R.normalize();
        
        double epsilon = 1e-3;
        Vector reflectedOrigin = intersection.position + epsilon * N;
        
        Ray reflectedRay(reflectedOrigin, R);
        
        Vector reflectedColor(0,0,0,0);
        double reflectedDepth = scene.camera.zFar;
        
        // APPEL RÉCURSIF : On lance un rayon dans la direction du reflet
        if (trace(reflectedRay, rayDepth, scene, reflectedColor, reflectedDepth)) {
            reflectedLight = reflectedColor;
        }
    }
    
    // [EXPL] RÉFRACTION (Transparence / Verre)
    Vector refractedLight(0);
    // refraction
    if ((material.refract > 1e-6) && (rayDepth < MAX_RAY_RECURSION))
    {
        // Loi de Snell-Descartes pour calculer le vecteur réfracté T
        // ... (calculs complexes des indices n1/n2 et du cosinus) ...
        Vector D = ray.direction.normalized();
        Vector N = intersection.normal;
        
        double cosI = D.dot(N);
        double n1 = 1.0; 
        double n2 = material.ior; 
        Vector normal = N;
        
        if (cosI > 0) {
            std::swap(n1, n2);
            normal = -N;
            cosI = -cosI;
        } else {
            cosI = -cosI;
        }
        
        double eta = n1 / n2;
        // total internal reflection check
        double k = 1.0 - eta * eta * (1.0 - cosI * cosI);
        
        if (k >= 0.0) { // Si pas de réflexion totale interne
            Vector T = eta * D + (eta * cosI - sqrt(k)) * normal;
            T.normalize();
            
            double epsilon = 1e-3;
            Vector refractedOrigin = intersection.position - epsilon * normal;
            
            Ray refractedRay(refractedOrigin, T);
            
            Vector refractedColor(0,0,0,0);
            double refractedDepth = scene.camera.zFar;
            // APPEL RÉCURSIF : On lance un rayon à travers l'objet
            if (trace(refractedRay, rayDepth, scene, refractedColor, refractedDepth)) {
                refractedLight = refractedColor;
            }
        }
    }
    
    double localWeight = 1.0; 
    // [EXPL] COMBINAISON FINALE
    // Couleur = Emission + (Ambiante+Diffuse+Spéculaire) + Reflet + Réfraction
    Vector finalColor = material.emission + 
           localWeight * (ambient + diffuse + specular) + 
           material.reflect * reflectedLight + 
           material.refract * refractedLight;

    return finalColor;
}