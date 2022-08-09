#pragma once
/**
 * @file rtx_global_illumination.cpp
 *
 * @brief various classes for calculating a ray tracing scene with global illumination
 *
 * @ingroup SundryLT
 *
 * @author Kevin Peivareh
 * Contact: kevin@peivareh.com
 *
 */

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <GL/glut.h>


/**
 * Implementation of a color in RGB color space
 *
 * to simplify the usage of RGB color in C++
 *
 */
class RGBSpace
{
   public:
      float R, G, B;

      void set(float r, float g, float b);
      std::string print();
      void mult(float c);
      void mult(RGBSpace c);
      void add(RGBSpace c);
      void sub(RGBSpace c);
      void clamp();
};

/**
 * Implementation of a position in 3D space
 *
 * to simplify the usage of positions in 3D space in C++
 *
 */
class Position
{
   public:
      float px, py, pz;

      void set(float x, float y, float z);
      std::string print();
      float distance(Position p);
};

/**
 * Implementation of a vector
 *
 * this class represents a vector which describes a direction in 3D space
 *
 */
class DirectionVector
{
   public:
      float vx, vy, vz;

      void set(float x, float y, float z);
      std::string print();
      void normalize();
      float dot(DirectionVector v);
      void mult(float c);
      void add(DirectionVector v);
      void sub(DirectionVector v);
};

/**
 * Implementation of light
 *
 * to simplify the usage of light rays in C++
 *
 */
class LightRay
{
   public:
      Position p;
      DirectionVector dir;

      void set(Position p, DirectionVector d);
      void set(Position p1, Position p2);
      std::string print();
      Position get_sample(float t);
};

/**
 * Implementation of a sphere in 3D space
 *
 * to simplify the usage of RGB color in C++
 *
 */
class Sphere
{
   public:
      Position center;
      DirectionVector motion;
      float radius;

      void set(Position p, float r) ;
      void set(Position p, DirectionVector m, float r) ;
      std::string print() ;
      bool calc_intersection(LightRay ray, Position &p, DirectionVector &normal);
};

/**
 * Implementation of PhongShading
 *
 * class for using the phong shading to shade different elements in 3D space
 *
 */
class PhongShading
{
 public:
   // Constructors
   PhongShading();
   ~PhongShading();

   // setter
   void InitCamera(Position pos);
   void LightColor(RGBSpace color, DirectionVector dir);
   void InitObject(RGBSpace color, float ka, float kd, float ks, float kp);
    
   // getter
   void GetShade(Position p, DirectionVector normal, RGBSpace & color);

 private:
   // cam
   Position CamPos;

   // light
   RGBSpace LightC;
   DirectionVector LightDir;
   
   // object
   RGBSpace ObjectColor;
   float Ka, Kd, Ks, Kp;
};