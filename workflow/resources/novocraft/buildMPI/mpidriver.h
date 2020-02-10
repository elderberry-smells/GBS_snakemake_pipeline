//
// File:   mpi_driver.h
// Author: Colin Hercus
// Copyright 2009-2015 Colin Hercus
//
// Created on 26 Nov 2009
//

#ifndef _mpi_driver_H
#define	_mpi_driver_H

#include <time.h>

/*
 *  MPI Message Buffer Size.
 * 
 * Larger Buffers  ...
 *    1. may increase MPI bandwidth and may mean you can use more nodes.
 *    2. may reduce concurrency at the end of the process hence causing a slow down.
 *    3. are more likely to help for fast alignments (small genomes, good quality reads, using --hlimit)
 *    4. are less likely to help for slow alignments (Biseq, Colour Space, 2 dye chemistry, not using --hlimit)
 */
static const int NOVO_MPI_MSG_LEN = 1024 * 512; // Try up to 1024 * 512;  
#ifdef MPITRACE
#define TRACE
#endif
#include "trace.h"
#define mpitracef tracef
#define mpitracel tracel
#define fdtracef tracef
#define fdtracel tracel
#define cstracef tracef
#define cstracel tracel

/**
 ** A generic class for MPI master and slave where only one type of work is being distributed to slaves
 ** All slave send results back to the master for collation.
 */
enum NV_Tag { TAG_NORMAL, TAG_EOF, TAG_MSGSIZE, TAG_SYNC, TAG_SYNCD, TAG_LAST };
enum NV_Datatype { DT_CHAR, DT_INT, DT_UNSIGNED, DT_SIZE };
typedef void * NV_MPI_STATUS;

struct mpiBuffer
{
    int size;
    int used;
    unsigned int tag;
    char message[0];
};
    
class mpi_driver {
public:
    static int rank;   // My MPI Rank or process number, 0 is the Master!
protected:
    int num_tasks; // The number of MPI tasks running.
    char *hostname;
    int l_hostname;
    int num_slave_threads;
    bool eof_received;

    // Dynamically allocated arrays for management of messages sent by master process.
    bool *eof_sent;
    bool *wait_sync;
    struct msgbuffers *masterRequest;
    // EOF Record, we should have our own unique datatype!
    static const unsigned char eof_message[4];
    static const int eof_len = 3;
    static const NV_Datatype eof_datatype;

    void run_slaves(int argc, char** argv);
    void run_master(int argc, char** argv);
public:
    mpi_driver();
    ~mpi_driver();
    static mpiBuffer * newMpiBuffer();
    static mpiBuffer * newMpiBuffer(int size);
    //Starts MPI Process and decides whether we are a master or slave
    void run_mpi(int argc, char** argv);
    // This is used to allocate buffers.
    int max_length();
    // Override with function to do any initialisation for the master process
    virtual void init_master(int argc, char** argv);
    // Override to construct a mesage to send;
    // Returns message length
    virtual mpiBuffer *next_message() = 0;
    void send_message(int task, int requestIdx);
    virtual void free_message(mpiBuffer **) = 0;
    
    virtual int master_process_message(int rank, char * buffer, int tag, int msg_len) = 0;
    // Returns true if no more messages to send.
    virtual bool eof() = 0;
    // Returns true if we need to sync the slaves.
    virtual bool syncPoint() { return false; };
    // Called when all slaves have ack'd the sync
    virtual void synced() { return; };

    // Called when all results have been received
    virtual bool final_results();

    // Create output processor
//    virtual mpi_collector* new_collector() = 0;

    // Slave routines
    virtual void init_slave(int argc, char** argv);

    virtual void process_message(int tag, int msg__len) = 0;
    // Called when eof is received.
    virtual void eof_slave();
    // Called after eof processed and barrier.
    virtual void final_slave();


};

int NV_MPI_Send(void* buf, int count, NV_Datatype datatype, int task, int tag);
int NV_MPI_Recv(void * stat_msg, int msg_len,  NV_Datatype datatype);
int NV_MPI_Recv(void * stat_msg, int msg_len,  NV_Datatype datatype, NV_MPI_STATUS);
int NV_MPI_Bcast(void* buf, int count, NV_Datatype datatype);
int NV_MPI_MsgLen(int count, NV_Datatype datatype);
NV_MPI_STATUS NV_Status_Object();
void deleteNV_Status_Object(NV_MPI_STATUS);
int NV_MPI_Get_count(NV_MPI_STATUS, NV_Datatype datatype);
int Barrier();
void Die();

#endif // _mpi_driver_H //
