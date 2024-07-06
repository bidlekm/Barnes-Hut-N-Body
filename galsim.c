#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>

double theta_max;
const double epsilon0 = 0.001;
const float circleRadius = 0.01, circleColor = 0;
const int windowWidth = 800;
double G = 100.0;

typedef struct Quad
{
    double topLeftX;
    double topLeftY;
    double botRightX;
    double botRightY;
    double centerOfMassX;
    double centerOfMassY;
    double *posX, *posY, *velX, *velY;
    double *particleMass;
    double quadMass;
    struct Quad *topLeftTree;
    struct Quad *topRightTree;
    struct Quad *botLeftTree;
    struct Quad *botRightTree;
    bool isLeaf;
} Quad;

void insert(Quad **quad, double *restrict posX, double *restrict posY, double *restrict velX, double *restrict velY,
            double *mass, const double topLeftX, const double topLeftY, const double botRightX, const double botRightY)
{
    if (*quad == NULL)
    {
        *quad = (Quad *)malloc(sizeof(Quad));
        (*quad)->topLeftX = topLeftX;
        (*quad)->topLeftY = topLeftY;
        (*quad)->botRightX = botRightX;
        (*quad)->botRightY = botRightY;

        (*quad)->posX = posX;
        (*quad)->posY = posY;
        (*quad)->velX = velX;
        (*quad)->velY = velY;
        (*quad)->particleMass = mass;
        (*quad)->isLeaf = true;

        (*quad)->topLeftTree = NULL;
        (*quad)->topRightTree = NULL;
        (*quad)->botLeftTree = NULL;
        (*quad)->botRightTree = NULL;
        return;
    }
    else
    {
        double midX = ((*quad)->topLeftX + (*quad)->botRightX) / 2;
        double midY = ((*quad)->topLeftY + (*quad)->botRightY) / 2;

        if ((*quad)->isLeaf == true)
        {
            // INSERT OLD PARTICLE
            double *oldParticleMass = (*quad)->particleMass;
            double *oldPosX = (*quad)->posX;
            double *oldPosY = (*quad)->posY;
            double *oldVelX = (*quad)->velX;
            double *oldVelY = (*quad)->velY;
            (*quad)->posX = NULL;
            (*quad)->posY = NULL;
            (*quad)->velX = NULL;
            (*quad)->velY = NULL;
            (*quad)->particleMass = NULL;

            if (*oldPosX <= midX)
            {
                if (*oldPosY <= midY)
                {
                    insert(&(*quad)->botLeftTree, oldPosX, oldPosY, oldVelX, oldVelY, oldParticleMass, (*quad)->topLeftX, midY, midX, (*quad)->botRightY);
                }
                else
                {
                    insert(&(*quad)->topLeftTree, oldPosX, oldPosY, oldVelX, oldVelY, oldParticleMass, (*quad)->topLeftX, (*quad)->topLeftY, midX, midY);
                }
            }
            else
            {
                if (*oldPosY <= midY)
                {
                    insert(&(*quad)->botRightTree, oldPosX, oldPosY, oldVelX, oldVelY, oldParticleMass, midX, midY, (*quad)->botRightX, (*quad)->botRightY);
                }
                else
                {
                    insert(&(*quad)->topRightTree, oldPosX, oldPosY, oldVelX, oldVelY, oldParticleMass, midX, (*quad)->topLeftY, (*quad)->botRightX, midY);
                }
            }
            (*quad)->isLeaf = false;
        }
        if (*posX <= midX)
        {
            if (*posY <= midY)
            {
                insert(&(*quad)->botLeftTree, posX, posY, velX, velY, mass, (*quad)->topLeftX, midY, midX, (*quad)->botRightY);
            }
            else
            {
                insert(&(*quad)->topLeftTree, posX, posY, velX, velY, mass, (*quad)->topLeftX, (*quad)->topLeftY, midX, midY);
            }
        }
        else
        {
            if (*posY <= midY)
            {
                insert(&(*quad)->botRightTree, posX, posY, velX, velY, mass, midX, midY, (*quad)->botRightX, (*quad)->botRightY);
            }
            else
            {
                insert(&(*quad)->topRightTree, posX, posY, velX, velY, mass, midX, (*quad)->topLeftY, (*quad)->botRightX, midY);
            }
        }
    }
}

double updateMass(Quad *restrict quad)
{
    // TODO: same thing as before, it was quad->particle
    if (quad->isLeaf)
    {
        quad->centerOfMassX = *quad->posX;
        quad->centerOfMassY = *quad->posY;
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

void calculateForcePair(const Quad *restrict quad, const double *restrict posX, const double *restrict posY, const double *restrict velX, const double *restrict velY, double *restrict forceX, double *restrict forceY, const double mass)
{
    double dx = *posX - quad->centerOfMassX;
    double dy = *posY - quad->centerOfMassY;
    double distSquare = (dx * dx) + (dy * dy);
    double dist = sqrt(distSquare) + epsilon0;
    double F_scalar = -G * mass * quad->quadMass / (dist * dist * dist);
    *forceX += F_scalar * dx;
    *forceY += F_scalar * dy;
}

void calculateForces(const Quad *restrict quad, const double *restrict posX, const double *restrict posY, const double *restrict velX, const double *restrict velY, double *restrict forceX, double *restrict forceY, const double mass)
{
    // TODO: SAME AS BEFORE
    if (quad->isLeaf && quad->posX != NULL && quad->posX != posX)
    {
        calculateForcePair(quad, posX, posY, velX, velY, forceX, forceY, mass);
    }
    else
    {
        double dx = *posX - quad->centerOfMassX;
        double dy = *posY - quad->centerOfMassY;
        double distSquare = (dx * dx) + (dy * dy);
        double theta = (quad->botRightX - quad->topLeftX) / sqrt(distSquare);
        if (theta <= theta_max)
        {
            calculateForcePair(quad, posX, posY, velX, velY, forceX, forceY, mass);
        }
        else
        {

            if (quad->topLeftTree != NULL)
            {
                calculateForces(quad->topLeftTree, posX, posY, velX, velY, forceX, forceY, mass);
            }
            if (quad->topRightTree != NULL)
            {
                calculateForces(quad->topRightTree, posX, posY, velX, velY, forceX, forceY, mass);
            }
            if (quad->botLeftTree != NULL)
            {
                calculateForces(quad->botLeftTree, posX, posY, velX, velY, forceX, forceY, mass);
            }
            if (quad->botRightTree != NULL)
            {
                calculateForces(quad->botRightTree, posX, posY, velX, velY, forceX, forceY, mass);
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

int readParticlesFromFile(const char *filename, double *restrict posX, double *restrict posY, double *restrict velX,
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

void printQuadtree(const Quad *quad, int depth)
{
    if (quad == NULL)
        return;

    if (quad->posX != NULL)
    {
        printf("Position: (%lf, %lf)\tMass: %lf\n", *quad->posX, *quad->posY, *(quad->particleMass));
    }
    if (quad->topLeftTree != NULL)
        printQuadtree(quad->topLeftTree, depth + 1);
    if (quad->topRightTree != NULL)
        printQuadtree(quad->topRightTree, depth + 1);
    if (quad->botLeftTree != NULL)
        printQuadtree(quad->botLeftTree, depth + 1);
    if (quad->botRightTree != NULL)
        printQuadtree(quad->botRightTree, depth + 1);
}

int main(int argc, char **argv)
{

    if (argc != 7)
    {
        printf("Not enough parameters, only %d provided\n", argc);
        return 1;
    }

    int N;
    char *filename = argv[2];
    int nsteps;
    double delta_t;
    int num_threads;
    if (sscanf(argv[1], "%d", &N) != 1 || sscanf(argv[3], "%d", &nsteps) != 1 ||
        sscanf(argv[4], "%lf", &delta_t) != 1 || sscanf(argv[5], "%lf", &theta_max) != 1 || sscanf(argv[6], "%d", &num_threads) != 1)
    {
        printf("Error with given parameters.\n");
        return 1;
    }
    omp_set_num_threads(num_threads);
    G = 100.0 / N;

    double *posX = (double *)malloc(N * sizeof(double));
    double *posY = (double *)malloc(N * sizeof(double));
    double *velX = (double *)malloc(N * sizeof(double));
    double *velY = (double *)malloc(N * sizeof(double));
    double *mass = (double *)malloc(N * sizeof(double));
    double *inv_mass = (double *)malloc(N * sizeof(double));
    double *brightness = (double *)malloc(N * sizeof(double));
    double *forcesX = (double *)calloc(N, sizeof(double));
    double *forcesY = (double *)calloc(N, sizeof(double));
    double *accelerationsX = (double *)calloc(N, sizeof(double));
    double *accelerationsY = (double *)calloc(N, sizeof(double));
    if (readParticlesFromFile(filename, posX, posY, velX, velY, mass, inv_mass, brightness, N) != 6 * N)
    {
        printf("Could not read in all the data from %s\n", filename);
        free(posX);
        free(posY);
        free(velX);
        free(velY);
        free(forcesX);
        free(forcesY);
        free(mass);
        free(inv_mass);
        free(accelerationsX);
        free(accelerationsY);
        free(brightness);
        return -1;
    }
    double startTime1 = omp_get_wtime();

    Quad *root = NULL;
    for (int i = 0; i < N; i++)
    {
        insert(&root, &posX[i], &posY[i], &velX[i], &velY[i], &mass[i], 0.0, 1.0, 1.0, 0.0);
    }
    updateMass(root);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; i++)
    {
        calculateForces(root, &posX[i], &posY[i], &velX[i], &velY[i], &forcesX[i], &forcesY[i], mass[i]);
    }
    destroyQuadtree(root);
    root = NULL;
    // MAIN LOOP
    for (int step = 0; step < nsteps; step++)
    {
        for (int i = 0; i < N; i++)
        {
            accelerationsX[i] = forcesX[i] * inv_mass[i];
            accelerationsY[i] = forcesY[i] * inv_mass[i];
            posX[i] += (delta_t * velX[i]) + (0.5 * delta_t * delta_t * accelerationsX[i]);
            posY[i] += (delta_t * velY[i]) + (0.5 * delta_t * delta_t * accelerationsY[i]);
        }
        memset(forcesX, 0, N * sizeof(double));
        memset(forcesY, 0, N * sizeof(double));
        Quad *root = NULL;

        for (int i = 0; i < N; i++)
        {
            insert(&root, &posX[i], &posY[i], &velX[i], &velY[i], &mass[i], 0.0, 1.0, 1.0, 0.0);
        }

        updateMass(root);

#pragma omp parallel for
        for (int i = 0; i < N; i++)
        {
            calculateForces(root, &posX[i], &posY[i], &velX[i], &velY[i], &forcesX[i], &forcesY[i], mass[i]);
        }

        for (int i = 0; i < N; i++)
        {
            velX[i] += 0.5 * delta_t * (forcesX[i] * inv_mass[i] + accelerationsX[i]);
            velY[i] += 0.5 * delta_t * (forcesY[i] * inv_mass[i] + accelerationsY[i]);
        }
        destroyQuadtree(root);
        root = NULL;
    }

    double totalTimeTaken = omp_get_wtime() - startTime1;

    writeParticlesToFile(posX, posY, velX, velY, mass, brightness, N);

    free(posX);
    free(posY);
    free(velX);
    free(velY);
    free(forcesX);
    free(forcesY);
    free(mass);
    free(inv_mass);
    free(accelerationsX);
    free(accelerationsY);
    free(brightness);

    printf("totalTimeTaken = %f\n", totalTimeTaken);
    return 0;
}