/**
 * @file TimeIntegrate.cpp
 * @author Luca Mazzotta (luca.mazzotta19@imperial.ac.uk)
 * @brief This file contains the function to perform the time integration.
 * @version 1.0
 * @date 2023-03-20
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "ShallowWater.h"
#include <omp.h>

using namespace std;

/**
 * @brief Performs the time integration using Runge-Kutta.
 * 
 */

void ShallowWater::TimeIntegrate()
    {
        // Initialising variables
        double *u = new double[Nx * Ny];
        double *v = new double[Nx * Ny];
        double *h = new double[Nx * Ny];

        double *dudx = new double[Nx * Ny];
        double *dudy = new double[Nx * Ny];

        double *dvdx = new double[Nx * Ny];
        double *dvdy = new double[Nx * Ny];

        double *dhdx = new double[Nx * Ny];
        double *dhdy = new double[Nx * Ny];

        double *k1u = new double[Nx * Ny];
        double *k1v = new double[Nx * Ny];
        double *k1h = new double[Nx * Ny];

        double *k2u = new double[Nx * Ny];
        double *k2v = new double[Nx * Ny];
        double *k2h = new double[Nx * Ny];

        double *k3u = new double[Nx * Ny];
        double *k3v = new double[Nx * Ny];
        double *k3h = new double[Nx * Ny];

        double *k4u = new double[Nx * Ny];
        double *k4v = new double[Nx * Ny];
        double *k4h = new double[Nx * Ny];

        double g = 9.81;

        for (int j = 0; j < T / dt; ++j)
        {
            // Copying U, V, H to u, v, h (temp arrays) and calculating hu, hv
            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                u[i] = U[i];
                v[i] = V[i];
                h[i] = H[i];
            }

            // Calculating K1
            xDerivative(u, dudx);
            yDerivative(u, dudy);
            xDerivative(v, dvdx);
            yDerivative(v, dvdy);
            xDerivative(h, dhdx);
            yDerivative(h, dhdy);

            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                k1u[i] = -u[i] * dudx[i] - v[i] * dudy[i] - g * dhdx[i];
                k1v[i] = -u[i] * dvdx[i] - v[i] * dvdy[i] - g * dhdy[i];
                k1h[i] = -(u[i] * dhdx[i] + h[i] * dudx[i]) - (v[i] * dhdy[i] + h[i] * dvdy[i]);
            }

            // Calculating K2
            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                u[i] = U[i] + k1u[i] * dt / 2;
                v[i] = V[i] + k1v[i] * dt / 2;
                h[i] = H[i] + k1h[i] * dt / 2;
            }

            xDerivative(u, dudx);
            yDerivative(u, dudy);
            xDerivative(v, dvdx);
            yDerivative(v, dvdy);
            xDerivative(h, dhdx);
            yDerivative(h, dhdy);

            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                k2u[i] = -u[i] * dudx[i] - v[i] * dudy[i] - g * dhdx[i];
                k2v[i] = -u[i] * dvdx[i] - v[i] * dvdy[i] - g * dhdy[i];
                k2h[i] = -(u[i] * dhdx[i] + h[i] * dudx[i]) - (v[i] * dhdy[i] + h[i] * dvdy[i]);
            }

            // Calculating K3
            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                u[i] = U[i] + k2u[i] * dt / 2;
                v[i] = V[i] + k2v[i] * dt / 2;
                h[i] = H[i] + k2h[i] * dt / 2;
            }

            xDerivative(u, dudx);
            yDerivative(u, dudy);
            xDerivative(v, dvdx);
            yDerivative(v, dvdy);
            xDerivative(h, dhdx);
            yDerivative(h, dhdy);

            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                k3u[i] = -u[i] * dudx[i] - v[i] * dudy[i] - g * dhdx[i];
                k3v[i] = -u[i] * dvdx[i] - v[i] * dvdy[i] - g * dhdy[i];
                k3h[i] = -(u[i] * dhdx[i] + h[i] * dudx[i]) - (v[i] * dhdy[i] + h[i] * dvdy[i]);
            }

            // Calculating K4
            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                u[i] = U[i] + k3u[i] * dt;
                v[i] = V[i] + k3v[i] * dt;
                h[i] = H[i] + k3h[i] * dt;
            }

            xDerivative(u, dudx);
            yDerivative(u, dudy);
            xDerivative(v, dvdx);
            yDerivative(v, dvdy);
            xDerivative(h, dhdx);
            yDerivative(h, dhdy);

            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                k4u[i] = -u[i] * dudx[i] - v[i] * dudy[i] - g * dhdx[i];
                k4v[i] = -u[i] * dvdx[i] - v[i] * dvdy[i] - g * dhdy[i];
                k4h[i] = -(u[i] * dhdx[i] + h[i] * dudx[i]) - (v[i] * dhdy[i] + h[i] * dvdy[i]);
            }

            // Calculating U, V, and H
            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                U[i] = U[i] + (k1u[i] + 2 * k2u[i] + 2 * k3u[i] + k4u[i]) * dt / 6;
                V[i] = V[i] + (k1v[i] + 2 * k2v[i] + 2 * k3v[i] + k4v[i]) * dt / 6;
                H[i] = H[i] + (k1h[i] + 2 * k2h[i] + 2 * k3h[i] + k4h[i]) * dt / 6;
            }
        }

        //Cleaning

        delete[] u;
        delete[] v;
        delete[] h;

        delete[] dudx;
        delete[] dudy;

        delete[] dvdx;
        delete[] dvdy;

        delete[] dhdx;
        delete[] dhdy;

        delete[] k1u;
        delete[] k1v;
        delete[] k1h;

        delete[] k2u;
        delete[] k2v;
        delete[] k2h;

        delete[] k3u;
        delete[] k3v;
        delete[] k3h;

        delete[] k4u;
        delete[] k4v;
        delete[] k4h;

    }