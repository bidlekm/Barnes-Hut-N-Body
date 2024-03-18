#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>

double theta_max = 0.5;
const double epsilon0 = 0.001;
const float circleRadius = 0.01, circleColor = 0;
const int windowWidth = 800;
double G = 100.0;

typedef struct Particle
{
    double *posX;
    double *posY;
    double *velX;
    double *velY;
} Particle;

typedef struct Quad
{
    double topLeftX;
    double topLeftY;
    double botRightX;
    double botRightY;
    double centerOfMassX;
    double centerOfMassY;
    Particle *particle;
    double *particleMass;
    double quadMass;
    struct Quad *topLeftTree;
    struct Quad *topRightTree;
    struct Quad *botLeftTree;
    struct Quad *botRightTree;
} Quad;

static double get_wall_seconds()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
}

void insert(Quad *quad, Particle *particle, double *mass)
{
    if (quad->particle == NULL && quad->topLeftTree == NULL && quad->topRightTree == NULL && quad->botLeftTree == NULL &&
        quad->botRightTree == NULL)
    {
        quad->particle = particle;
        quad->particleMass = mass;
        return;
    }
    else
    {
        double midX = (quad->topLeftX + quad->botRightX) / 2;
        double midY = (quad->topLeftY + quad->botRightY) / 2;

        if (*particle->posX <= midX)
        {
            if (*particle->posY <= midY)
            {
                if (quad->botLeftTree == NULL)
                {
                    Quad *newQuad = (Quad *)malloc(sizeof(Quad));
                    newQuad->topLeftX = quad->topLeftX;
                    newQuad->topLeftY = midY;
                    newQuad->botRightX = midX;
                    newQuad->botRightY = quad->botRightY;
                    newQuad->quadMass = 0;
                    newQuad->particle = NULL;
                    newQuad->particleMass = NULL;
                    newQuad->centerOfMassX = 0;
                    newQuad->centerOfMassY = 0;
                    newQuad->topLeftTree = NULL;
                    newQuad->botRightTree = NULL;
                    newQuad->topRightTree = NULL;
                    newQuad->botLeftTree = NULL;
                    quad->botLeftTree = newQuad;
                }
                insert(quad->botLeftTree, particle, mass);
            }
            else
            {
                if (quad->topLeftTree == NULL)
                {
                    Quad *newQuad = (Quad *)malloc(sizeof(Quad));
                    newQuad->topLeftX = quad->topLeftX;
                    newQuad->topLeftY = quad->topLeftY;
                    newQuad->botRightX = midX;
                    newQuad->botRightY = midY;
                    newQuad->quadMass = 0;
                    newQuad->particle = NULL;
                    newQuad->particleMass = NULL;
                    newQuad->centerOfMassX = 0;
                    newQuad->centerOfMassY = 0;
                    newQuad->topLeftTree = NULL;
                    newQuad->botRightTree = NULL;
                    newQuad->topRightTree = NULL;
                    newQuad->botLeftTree = NULL;
                    quad->topLeftTree = newQuad;
                }
                insert(quad->topLeftTree, particle, mass);
            }
        }
        else
        {
            if (*particle->posY <= midY)
            {
                if (quad->botRightTree == NULL)
                {
                    Quad *newQuad = (Quad *)malloc(sizeof(Quad));
                    newQuad->topLeftX = midX;
                    newQuad->topLeftY = midY;
                    newQuad->botRightX = quad->botRightX;
                    newQuad->botRightY = quad->botRightY;
                    newQuad->quadMass = 0;
                    newQuad->particle = NULL;
                    newQuad->particleMass = NULL;
                    newQuad->centerOfMassX = 0;
                    newQuad->centerOfMassY = 0;
                    newQuad->topLeftTree = NULL;
                    newQuad->botRightTree = NULL;
                    newQuad->topRightTree = NULL;
                    newQuad->botLeftTree = NULL;
                    quad->botRightTree = newQuad;
                }
                insert(quad->botRightTree, particle, mass);
            }
            else
            {
                if (quad->topRightTree == NULL)
                {
                    Quad *newQuad = (Quad *)malloc(sizeof(Quad));
                    newQuad->topLeftX = midX;
                    newQuad->topLeftY = quad->topLeftY;
                    newQuad->botRightX = quad->botRightX;
                    newQuad->botRightY = midY;
                    newQuad->quadMass = 0;
                    newQuad->particle = NULL;
                    newQuad->particleMass = NULL;
                    newQuad->centerOfMassX = 0;
                    newQuad->centerOfMassY = 0;
                    newQuad->topLeftTree = NULL;
                    newQuad->botRightTree = NULL;
                    newQuad->topRightTree = NULL;
                    newQuad->botLeftTree = NULL;
                    quad->topRightTree = newQuad;
                }
                insert(quad->topRightTree, particle, mass);
            }
        }

        if (quad->particle != NULL)
        {
            Particle *oldParticle = quad->particle;
            double *oldParticleMass = quad->particleMass;
            quad->particle = NULL;
            quad->particleMass = NULL;
            insert(quad, oldParticle, oldParticleMass);
        }
    }
}

double updateMass(Quad *restrict quad)
{
    if (quad->particle != NULL)
    {
        quad->centerOfMassX = *quad->particle->posX;
        quad->centerOfMassY = *quad->particle->posY;
        quad->quadMass = *quad->particleMass;
        return quad->quadMass;
    }

    double totalMass = 0;
    double centerX = 0;
    double centerY = 0;

    if (quad->topLeftTree != NULL)
    {
        double mass = updateMass(quad->topLeftTree);
        centerX += mass * quad->topLeftTree->centerOfMassX;
        centerY += mass * quad->topLeftTree->centerOfMassY;
        totalMass += mass;
    }
    if (quad->topRightTree != NULL)
    {
        double mass = updateMass(quad->topRightTree);
        centerX += mass * quad->topRightTree->centerOfMassX;
        centerY += mass * quad->topRightTree->centerOfMassY;
        totalMass += mass;
    }
    if (quad->botLeftTree != NULL)
    {
        double mass = updateMass(quad->botLeftTree);
        centerX += mass * quad->botLeftTree->centerOfMassX;
        centerY += mass * quad->botLeftTree->centerOfMassY;
        totalMass += mass;
    }
    if (quad->botRightTree != NULL)
    {
        double mass = updateMass(quad->botRightTree);
        centerX += mass * quad->botRightTree->centerOfMassX;
        centerY += mass * quad->botRightTree->centerOfMassY;
        totalMass += mass;
    }

    quad->quadMass = totalMass;
    quad->centerOfMassX = centerX / totalMass;
    quad->centerOfMassY = centerY / totalMass;

    return totalMass;
}

void calculateForcePair(const Quad *restrict quad, const Particle *restrict particle, double *restrict forceX, double *restrict forceY, const double mass)
{
    double dx = *particle->posX - quad->centerOfMassX;
    double dy = *particle->posY - quad->centerOfMassY;
    double distSquare = (dx * dx) + (dy * dy);
    double dist = sqrt(distSquare) + epsilon0;
    double F_scalar = -G * mass * quad->quadMass / (dist * dist * dist);
    *forceX += F_scalar * dx;
    *forceY += F_scalar * dy;
}

void calculateForces(const Quad *restrict quad, const Particle *restrict particle, double *restrict forceX, double *restrict forceY, const double mass)
{
    if (quad->topLeftTree == NULL && quad->topRightTree == NULL && quad->botLeftTree == NULL &&
        quad->botRightTree == NULL && quad->particle != NULL && quad->particle != particle)
    {
        calculateForcePair(quad, particle, forceX, forceY, mass);
    }
    else
    {
        double dx = *particle->posX - quad->centerOfMassX;
        double dy = *particle->posY - quad->centerOfMassY;
        double distSquare = (dx * dx) + (dy * dy);
        double theta = (quad->botRightX - quad->topLeftX) / sqrt(distSquare);
        if (theta <= theta_max)
        {
            calculateForcePair(quad, particle, forceX, forceY, mass);
        }
        else
        {
            if (quad->topLeftTree != NULL)
            {
                calculateForces(quad->topLeftTree, particle, forceX, forceY, mass);
            }
            if (quad->topRightTree != NULL)
            {
                calculateForces(quad->topRightTree, particle, forceX, forceY, mass);
            }
            if (quad->botLeftTree != NULL)
            {
                calculateForces(quad->botLeftTree, particle, forceX, forceY, mass);
            }
            if (quad->botRightTree != NULL)
            {
                calculateForces(quad->botRightTree, particle, forceX, forceY, mass);
            }
        }
    }
}

void destroyQuadtree(Quad *quad)
{
    if (quad == NULL)
    {
        return;
    }
    destroyQuadtree(quad->topLeftTree);
    destroyQuadtree(quad->topRightTree);
    destroyQuadtree(quad->botLeftTree);
    destroyQuadtree(quad->botRightTree);
    free(quad);
    quad = NULL;
}

int readParticlesFromFile(const char *filename, Particle *particles, double *restrict posX, double *restrict posY, double *restrict velX,
                          double *restrict velY, double *restrict mass, double *restrict inv_mass, double *restrict brightness, const int N)
{
    FILE *inputFile = fopen(filename, "rb");
    if (inputFile == NULL)
    {
        printf("%s not found\n", filename);
        return 0;
    }
    int reads = 0;
    for (int i = 0; i < N; i++)
    {
        reads += fread(&(posX[i]), sizeof(double), 1, inputFile);
        reads += fread(&(posY[i]), sizeof(double), 1, inputFile);
        reads += fread(&(mass[i]), sizeof(double), 1, inputFile);
        reads += fread(&(velX[i]), sizeof(double), 1, inputFile);
        reads += fread(&(velY[i]), sizeof(double), 1, inputFile);
        reads += fread(&(brightness[i]), sizeof(double), 1, inputFile);
        particles[i].posX = &posX[i];
        particles[i].posY = &posY[i];
        particles[i].velX = &velX[i];
        particles[i].velY = &velY[i];

        inv_mass[i] = 1 / mass[i];
    }
    fclose(inputFile);
    return reads;
}

void writeParticlesToFile(double *restrict posX, double *restrict posY, double *restrict velX, double *restrict velY, const double *restrict mass,
                          const double *restrict brightness, const int N)
{
    FILE *outputFile = fopen("result.gal", "wb");
    if (outputFile == NULL)
    {
        perror("Error opening file");
        return;
    }

    for (int i = 0; i < N; i++)
    {
        fwrite(&(posX[i]), sizeof(double), 1, outputFile);
        fwrite(&(posY[i]), sizeof(double), 1, outputFile);
        fwrite(&(mass[i]), sizeof(double), 1, outputFile);
        fwrite(&(velX[i]), sizeof(double), 1, outputFile);
        fwrite(&(velY[i]), sizeof(double), 1, outputFile);
        fwrite(&(brightness[i]), sizeof(double), 1, outputFile);
    }

    fclose(outputFile);
    // printf("Particles written to result.gal successfully\n");
}

void printAllParticles(double *restrict posX, double *restrict posY, double *restrict velX, double *restrict velY,
                       const double *mass, const double *brightness, const int N)
{
    for (int i = 0; i < N; i++)
    {
        printf("Particle %d: position %lf, %lf; mass %lf; velocity %lf, %lf; brightness %lf\n", i, posX[i], posY[i],
               mass[i], velX[i], velY[i], brightness[i]);
    }
}

int main(int argc, char **argv)
{
    double startTime1 = get_wall_seconds();
    if (argc != 5)
    {
        printf("Not enough parameters, only %d provided\n", argc);
        return 1;
    }

    int N;
    char *filename = argv[2];
    int nsteps;
    double delta_t;
    if (sscanf(argv[1], "%d", &N) != 1 || sscanf(argv[3], "%d", &nsteps) != 1 ||
        sscanf(argv[4], "%lf", &delta_t) != 1)
    {
        printf("Error with given parameters.\n");
        return 1;
    }
    omp_set_num_threads(8);
    G = 100.0 / N;

    Particle *particles = (Particle *)malloc(N * sizeof(Particle));
    double *posX = (double *)malloc(N * sizeof(Particle));
    double *posY = (double *)malloc(N * sizeof(Particle));
    double *velX = (double *)malloc(N * sizeof(Particle));
    double *velY = (double *)malloc(N * sizeof(Particle));
    double *mass = (double *)malloc(N * sizeof(double));
    double *inv_mass = (double *)malloc(N * sizeof(double));
    double *brightness = (double *)malloc(N * sizeof(double));
    double *forcesX = (double *)calloc(N, sizeof(double));
    double *forcesY = (double *)calloc(N, sizeof(double));
    double *accelerationsX = (double *)calloc(N, sizeof(double));
    double *accelerationsY = (double *)calloc(N, sizeof(double));
    if (readParticlesFromFile(filename, particles, posX, posY, velX, velY, mass, inv_mass, brightness, N) != 6 * N)
    {
        printf("Could not read in all the data from %s\n", filename);
        free(posX);
        free(posY);
        free(velX);
        free(velY);
        free(particles);
        free(forcesX);
        free(forcesY);
        free(mass);
        free(inv_mass);
        free(accelerationsX);
        free(accelerationsY);
        free(brightness);
        return -1;
    }

    Quad *root = (Quad *)malloc(sizeof(Quad));
    root->topLeftX = 0.0;
    root->topLeftY = 1.0;
    root->botRightX = 1.0;
    root->botRightY = 0.0;
    root->quadMass = 0;
    root->centerOfMassX = 0;
    root->centerOfMassY = 0;
    root->particle = NULL;
    root->particleMass = NULL;
    root->topLeftTree = NULL;
    root->botRightTree = NULL;
    root->botLeftTree = NULL;
    root->topRightTree = NULL;
    for (int i = 0; i < N; i++)
    {
        insert(root, &particles[i], &mass[i]);
    }
    updateMass(root);
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        calculateForces(root, &particles[i], &forcesX[i], &forcesY[i], mass[i]);
    }
    destroyQuadtree(root);

    for (int step = 0; step < nsteps; step++)
    {
#pragma omp parallel for
        for (int i = 0; i < N; i++)
        {
            accelerationsX[i] = forcesX[i] * inv_mass[i];
            accelerationsY[i] = forcesY[i] * inv_mass[i];
            *particles[i].posX += delta_t * (*particles[i].velX) + 0.5 * delta_t * delta_t * accelerationsX[i];
            *particles[i].posY += delta_t * (*particles[i].velY) + 0.5 * delta_t * delta_t * accelerationsY[i];
        }
        memset(forcesX, 0, N * sizeof(double));
        memset(forcesY, 0, N * sizeof(double));
        Quad *root = (Quad *)malloc(sizeof(Quad));
        root->topLeftX = 0.0;
        root->topLeftY = 1.0;
        root->botRightX = 1.0;
        root->botRightY = 0.0;
        root->quadMass = 0;
        root->centerOfMassX = 0;
        root->centerOfMassY = 0;
        root->particle = NULL;
        root->particleMass = NULL;
        root->topLeftTree = NULL;
        root->botRightTree = NULL;
        root->botLeftTree = NULL;
        root->topRightTree = NULL;

        for (int i = 0; i < N; i++)
        {
            insert(root, &particles[i], &mass[i]);
        }

        updateMass(root);

#pragma omp parallel for
        for (int i = 0; i < N; i++)
        {
            calculateForces(root, &particles[i], &forcesX[i], &forcesY[i], mass[i]);
        }
#pragma omp parallel for
        for (int i = 0; i < N; i++)
        {
            *particles[i].velX += 0.5 * delta_t * (forcesX[i] * inv_mass[i] + accelerationsX[i]);
            *particles[i].velY += 0.5 * delta_t * (forcesY[i] * inv_mass[i] + accelerationsY[i]);
        }

        destroyQuadtree(root);
        root = NULL;
    }

    writeParticlesToFile(posX, posY, velX, velY, mass, brightness, N);

    free(posX);
    free(posY);
    free(velX);
    free(velY);
    free(particles);
    free(forcesX);
    free(forcesY);
    free(mass);
    free(inv_mass);
    free(accelerationsX);
    free(accelerationsY);
    free(brightness);

    double totalTimeTaken = get_wall_seconds() - startTime1;
    printf("totalTimeTaken = %f\n", totalTimeTaken);
    return 0;
}
