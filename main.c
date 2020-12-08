#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

typedef struct NeighborRank {
    int up;
    int down;
    int left;
    int right;
} NeighborRank;

#define MASTER_RANK     0                  
#define GENERATIONS  1000                

#define reduce_enabled()    (REDUCE != 0)
#define is_master() (rank == MASTER_RANK)

static int alive = 0;
static int changed = 0;

/*---------------------- Inline sunartiseis ---------------------------*/

inline int updateInternal(int columns, int rows, char *u1, char *u2);
inline int updateBorders(int columns, int rows, char *u1, char *u2);
inline int calculateNextState(int columns, int rows, char *u1, char *u2, int x, int y);

/*---------------------- Main ---------------------------*/

int life(const int rank, const int N, MPI_Comm * cartComm, const int GRID_SIDE, const int REDUCE) {
    MPI_Pcontrol(0);

    NeighborRank neighbors = {0};
    const int TOPOLOGY_GRID_SIDE_HEIGHT = (int) (sqrt(N));
    const int TOPOLOGY_GRID_SIDE_WIDTH = N / TOPOLOGY_GRID_SIDE_HEIGHT;
    const int TOPOLOGY_DIMENSIONS[2] = {TOPOLOGY_GRID_SIDE_HEIGHT, TOPOLOGY_GRID_SIDE_WIDTH};
    const int TOPOLOGY_PERIOD[2] = {0, 0};
    const int PROCESS_COLUMNS = GRID_SIDE / TOPOLOGY_GRID_SIDE_WIDTH;
    const int PROCESS_ROWS = GRID_SIDE / TOPOLOGY_GRID_SIDE_HEIGHT;
    const int HALO_COLUMNS = PROCESS_COLUMNS + 2;
    const int HALO_ROWS = PROCESS_ROWS + 2;
    const int GRID_CELLS = GRID_SIDE*GRID_SIDE;
    const int HALO_CELLS = HALO_COLUMNS * HALO_ROWS;
    double start_time, end_time;

    MPI_Request incomingMessage[2][4], outgoingMessage[2][4];

   


    printf("Process R-%d, block rows: %d block cols:%d \n", rank, PROCESS_ROWS, PROCESS_COLUMNS);

    if (MPI_Cart_create(MPI_COMM_WORLD, 2, TOPOLOGY_DIMENSIONS, TOPOLOGY_PERIOD, 1, cartComm) != MPI_SUCCESS) {
        printf("MPI_SUCCESS not returned in MPI_Cart_create ! \n");
        exit(2);
    } else {
        printf("MPI_CART created with dimensions: %dx%d \n", TOPOLOGY_GRID_SIDE_HEIGHT, TOPOLOGY_GRID_SIDE_WIDTH);
    }

    if (MPI_Cart_shift(*cartComm, 0, 1, &neighbors.up, &neighbors.down) != MPI_SUCCESS) {
        printf("MPI_SUCCESS not returned in MPI_Cart_shift ! \n");
    }
    if (MPI_Cart_shift(*cartComm, 1, 1, &neighbors.left, &neighbors.right) != MPI_SUCCESS) {
        printf("MPI_SUCCESS not returned in MPI_Cart_shift ! \n");
    }

    printf("Hello from process with rank: %d, neighbors: %d %d %d %d \n", rank, neighbors.up, neighbors.right, neighbors.down, neighbors.left);

    /*-------------------------------------- datatypes -------------------------------------*/

    MPI_Datatype workerInputDatatype, workerOutputDatatype, columnDatatype, rowDatatype;

    // custom datatypes gia stiles kai seires

    if (MPI_Type_vector(PROCESS_ROWS, 1, HALO_COLUMNS, MPI_BYTE, &columnDatatype) != MPI_SUCCESS) {
        printf("MPI_SUCCESS not returned in MPI_Type_vector ! \n");
        exit(1);
    }

    if (MPI_Type_vector(1, PROCESS_COLUMNS, HALO_COLUMNS, MPI_BYTE, &rowDatatype) != MPI_SUCCESS) {
        printf("MPI_SUCCESS not returned in MPI_Type_vector ! \n");
        exit(1);
    }

    if (MPI_Type_vector(PROCESS_ROWS, PROCESS_COLUMNS, GRID_SIDE, MPI_BYTE, &workerInputDatatype) != MPI_SUCCESS) {
        printf("MPI_SUCCESS not returned in MPI_Type_vector ! \n");
        exit(1);
    }

    if (MPI_Type_vector(PROCESS_ROWS, PROCESS_COLUMNS, HALO_COLUMNS, MPI_BYTE, &workerOutputDatatype) != MPI_SUCCESS) {
        printf("MPI_SUCCESS not returned in MPI_Type_vector ! \n");
        exit(1);
    }

    if (MPI_Type_commit(&columnDatatype) != MPI_SUCCESS) {
        printf("MPI_SUCCESS not returned in MPI_Type_commit ! \n");
        exit(1);
    }

    if (MPI_Type_commit(&rowDatatype) != MPI_SUCCESS) {
        printf("MPI_SUCCESS not returned in MPI_Type_commit ! \n");
        exit(1);
    }

    if (MPI_Type_commit(&workerInputDatatype) != MPI_SUCCESS) {
        printf("MPI_SUCCESS not returned in masterToWorker ! \n");
        exit(1);
    }

    if (MPI_Type_commit(&workerOutputDatatype) != MPI_SUCCESS) {
        printf("MPI_SUCCESS not returned in workerToMaster ! \n");
        exit(1);
    }


    char * currentGen = malloc(HALO_CELLS * sizeof (char));
    char * nextGen = malloc(HALO_CELLS * sizeof (char));

    for (int i = 0; i < HALO_CELLS; i++) {
        currentGen[i] = rand() % 2;
    }

    // ------------------------------------------------------------------
    int offsetX, offsetY; 


    int swap; 


    char * temp = currentGen;
    char temp_i = 0;

    if (MPI_Send_init(&temp[HALO_COLUMNS + 1], 1, rowDatatype, neighbors.up, 0, *cartComm, &outgoingMessage[temp_i][0]) != MPI_SUCCESS) {
        printf("Error in creating request Send 0 - 0 \n");
    }

    if (MPI_Send_init(&temp[2 * HALO_COLUMNS - 2], 1, columnDatatype, neighbors.right, 0, *cartComm, &outgoingMessage[temp_i][2]) != MPI_SUCCESS) {
        printf("Error in creating request Send 0 - 2 \n");
    }
    if (MPI_Send_init(&temp[(HALO_ROWS - 2) * HALO_COLUMNS + 1], 1, rowDatatype, neighbors.down, 0, *cartComm, &outgoingMessage[temp_i][1]) != MPI_SUCCESS) {
        printf("Error in creating request Send 0 - 1 \n");
    }
    if (MPI_Send_init(&temp[HALO_COLUMNS + 1], 1, columnDatatype, neighbors.left, 0, *cartComm, &outgoingMessage[temp_i][3]) != MPI_SUCCESS) {
        printf("Error in creating request Send 0 - 3 \n");
    }

    if (MPI_Recv_init(&temp[1], 1, rowDatatype, neighbors.up, 0, *cartComm, &incomingMessage[temp_i][0]) != MPI_SUCCESS) {
        printf("Error in creating request Recv 0 - 0 \n");
    }
    if (MPI_Recv_init(&temp[2 * HALO_COLUMNS - 1], 1, columnDatatype, neighbors.right, 0, *cartComm, &incomingMessage[temp_i][2]) != MPI_SUCCESS) {
        printf("Error in creating request Recv 0 - 2 \n");
    }
    if (MPI_Recv_init(&temp[(HALO_ROWS - 1) * HALO_COLUMNS + 1], 1, rowDatatype, neighbors.down, 0, *cartComm, &incomingMessage[temp_i][1]) != MPI_SUCCESS) {
        printf("Error in creating request Recv 0 - 1 \n");
    }
    if (MPI_Recv_init(&temp[HALO_COLUMNS], 1, columnDatatype, neighbors.left, 0, *cartComm, &incomingMessage[temp_i][3]) != MPI_SUCCESS) {
        printf("Error in creating request Recv 0 - 3 \n");
    }

    temp = nextGen;
    temp_i = 1;

    if (MPI_Send_init(&temp[HALO_COLUMNS + 1], 1, rowDatatype, neighbors.up, 0, *cartComm, &outgoingMessage[temp_i][0]) != MPI_SUCCESS) {
        printf("Error in creating request Send 1 - 0 \n");
    }
    if (MPI_Send_init(&temp[2 * HALO_COLUMNS - 2], 1, columnDatatype, neighbors.right, 0, *cartComm, &outgoingMessage[temp_i][2]) != MPI_SUCCESS) {
        printf("Error in creating request Send 1 - 0 \n");
    }
    if (MPI_Send_init(&temp[(HALO_ROWS - 2) * HALO_COLUMNS + 1], 1, rowDatatype, neighbors.down, 0, *cartComm, &outgoingMessage[temp_i][1]) != MPI_SUCCESS) {
        printf("Error in creating request Send 1 - 0 \n");
    }
    if (MPI_Send_init(&temp[HALO_COLUMNS + 1], 1, columnDatatype, neighbors.left, 0, *cartComm, &outgoingMessage[temp_i][3]) != MPI_SUCCESS) {
        printf("Error in creating request Send 1 - 0 \n");
    }

    if (MPI_Recv_init(&temp[1], 1, rowDatatype, neighbors.up, 0, *cartComm, &incomingMessage[temp_i][0]) != MPI_SUCCESS) {
        printf("Error in creating request Recv 1 - 0 \n");
    }
    if (MPI_Recv_init(&temp[2 * HALO_COLUMNS - 1], 1, columnDatatype, neighbors.right, 0, *cartComm, &incomingMessage[temp_i][2]) != MPI_SUCCESS) {
        printf("Error in creating request Recv 1 - 2 \n");
    }
    if (MPI_Recv_init(&temp[(HALO_ROWS - 1) * HALO_COLUMNS + 1], 1, rowDatatype, neighbors.down, 0, *cartComm, &incomingMessage[temp_i][1]) != MPI_SUCCESS) {
        printf("Error in creating request Recv 1 - 1 \n");
    }
    if (MPI_Recv_init(&temp[HALO_COLUMNS], 1, columnDatatype, neighbors.left, 0, *cartComm, &incomingMessage[temp_i][3]) != MPI_SUCCESS) {
        printf("Error in creating request Recv 1 - 3 \n");
    }

    swap = 0;
    MPI_Barrier(*cartComm);

    int reduce1, reduce2, isReduce, reduction;

    printf("Rank %d received work. Beginning time steps...\n", rank);

    if (rank == MASTER_RANK) {
        start_time = MPI_Wtime();
    }


    MPI_Pcontrol(1);

    for (int it = 1; it <= GENERATIONS; it++) {
        alive = 0;
        changed = 0;

        if (MPI_Startall(4, &incomingMessage[swap][0]) != MPI_SUCCESS) {
            printf("MPI_Startall failed on incoming message: %d %d %d \n", it, swap, 0);
        }

        if (MPI_Startall(4, &outgoingMessage[swap][0]) != MPI_SUCCESS) {
            printf("MPI_Startall failed on outgoing message: %d %d %d \n", it, swap, 0);
        }

        if (updateInternal(HALO_COLUMNS, HALO_ROWS, currentGen, nextGen) != MPI_SUCCESS) {
            printf("updateInternal failed on outgoing message: %d %d %d \n", it, swap, 0);
        }

        if (MPI_Waitall(4, &incomingMessage[swap][0], MPI_STATUSES_IGNORE) != MPI_SUCCESS) {
            printf("MPI_Waitall failed on incoming message: %d %d %d \n", it, swap, 0);
        }

        if (updateBorders(HALO_COLUMNS, HALO_ROWS, currentGen, nextGen) != MPI_SUCCESS) {
            printf("updateBorders failed on outgoing message: %d %d %d \n", it, swap, 0);
        }

        if (MPI_Waitall(4, &outgoingMessage[swap][0], MPI_STATUSES_IGNORE) != MPI_SUCCESS) {
            printf("MPI_Waitall failed on incoming message: %d %d %d \n", it, swap, 0);
        }


        temp = currentGen; // Swapping pinakwn anti gia copy
        currentGen = nextGen;
        nextGen = temp;
        swap = 1 - swap;

        //reduction

        if (reduce_enabled() && it % REDUCE == 0) {
            isReduce = alive == 0 || changed == 0; 
            MPI_Allreduce(&isReduce, &reduction, 1, MPI_INT, MPI_LOR, *cartComm); 
        }
    }

    MPI_Pcontrol(0);

    

    if (rank == MASTER_RANK) { 
        end_time = MPI_Wtime();
        printf("Elapsed time = %f seconds\n", end_time - start_time); //print total computation time
    }

    free(currentGen);
    free(nextGen);

    return 0;
}

int main(int argc, char *argv[]) {
    MPI_Comm cartComm = {0};

    int REDUCE = 0;
    int GRID_SIDE = 420;
    int N = 0;
    int rank = 0; // Process/Task unique id

    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--reduce") == 0 || strcmp(argv[i], "-r") == 0) {
            REDUCE = 1;
        }

        if (strcmp(argv[i], "--grid_size") == 0 || strcmp(argv[i], "-g") == 0) {
            GRID_SIDE = atoi(argv[i + 1]);
        }
    }
    
    MPI_Init(&argc, &argv); // Initialize mpi
    MPI_Comm_size(MPI_COMM_WORLD, &N); // anathesi tasks stin topologia
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // anathesi ranks stin topologia

    if (is_master()) {
        const int SEED = time(NULL);
        srand(SEED);

        printf("Random seed             : %d \n", SEED);
        printf("Grid single dimension   : %d \n", GRID_SIDE);
        printf("Grid dimensions         : %d x %d \n", GRID_SIDE, GRID_SIDE);
        printf("Halo dimensions         : %d x %d \n", GRID_SIDE + 2, GRID_SIDE + 2);

        if (REDUCE != 0) {
            printf("Reduce mode            : %d \n", REDUCE);
        } else {
            printf("Reduce mode            : off \n");
        }
    }

    life(rank, N, &cartComm, GRID_SIDE, REDUCE);

    MPI_Finalize();

    return 0;
}

int calculateNextState(int columns, int rows, char *u1, char *u2, int x, int y) {
    int coordinate1D = y * columns + x;

    int count = -u1[coordinate1D];
    for (int i = -1; i <= +1; i++) { // -1, 0, 1
        for (int j = -1; j <= +1; j++) { // -1 ,0 ,1
            count = count + u1[coordinate1D + (i * columns) + j];
        }
    }
    //printf(" y:%d x:%d count:%d \n",y,x,count);
    if (count < 2 || count > 3) {
        u2[coordinate1D] = 0;
    } else if (count == 2) {
        u2[coordinate1D] = u1[coordinate1D];
        if (u2[coordinate1D] == 1) {
            alive = 1;
        }
    } else {
        u2[coordinate1D] = 1;
        alive = 1;
    }

    if (u2[coordinate1D] != u1[coordinate1D]) {
        changed = 1;
    }
}

/*----------------------- sunartiseis enimerwsis --------------------------*/

int updateInternal(int columns, int rows, char *u1, char *u2) {
    for (int y = 2; y < rows - 2; y++) {
        for (int x = 2; x < columns - 2; x++) {
            calculateNextState(columns, rows, u1, u2, x, y);
        }
    }

    return MPI_SUCCESS;
}

int updateBorders(int columns, int rows, char *u1, char *u2) {
    int x, y;

    for (int x = 1, y = 1; x <= columns - 2; x++) {
        calculateNextState(columns, rows, u1, u2, x, y);
    }

    for (int x = 1, y = rows - 2; x <= columns - 2; x++) {
        calculateNextState(columns, rows, u1, u2, x, y);
    }

    for (int y = 1, x = 1; y <= rows - 2; y++) {
        calculateNextState(columns, rows, u1, u2, x, y);
    }

    for (int y = 1, x = columns - 2; y <= rows - 2; y++) {
        calculateNextState(columns, rows, u1, u2, x, y);
    }

    return MPI_SUCCESS;
}
