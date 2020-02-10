//
// File:   mpi_driver.cc
// Author: Colin Hercus
// Copyright 2009 Colin Hercus
//
// Created on 26 Nov 2009
//

#include <stdio.h>
#include <stdlib.h>
#include "mpidriver.h"
#include <mpi.h>
#include "TTimer.hh"
#include "novolib.h"

const unsigned char mpi_driver::eof_message[4] = "EOF";
const NV_Datatype mpi_driver::eof_datatype = DT_CHAR;
MPI_Datatype dtConvert[DT_SIZE] = {MPI_CHAR, MPI_INT, MPI_UNSIGNED};

int mpi_driver::rank = 0;

struct msgbuffers {
    int size;
    MPI_Request *requests;
    struct mpiBuffer **buffers;
    
    msgbuffers()
    {
        requests = NULL;
        buffers = NULL;
        size = 0;
    }

    ~msgbuffers()
    {
        mpitracel;
        delete[] requests;
        for(int i = 0; i < size; i++) 
        {
            mpitracef("%d", i);
            if(buffers[i] != NULL)
                free(buffers[i]);
        }
        delete[] buffers;
        mpitracef("Exit");
    }

    void allocate(int nRcv, int nTot)
    {
        mpitracef("%d:%d buffers of %d bytes", nRcv, nTot, NOVO_MPI_MSG_LEN);
        size = nTot;
        requests = new MPI_Request[size];
        buffers = new struct mpiBuffer *[size];
        int i = 0;
        for(; i < nRcv; i++)
        {
            requests[i] = MPI_REQUEST_NULL;
            buffers[i] = mpi_driver::newMpiBuffer(); 
        }
        for(; i < nTot; i++)
        {
            requests[i] = MPI_REQUEST_NULL;
            buffers[i] = NULL; 
        }
    }

    void resize(int i, unsigned int bufsize) {
        if(buffers[i] != NULL && bufsize < buffers[i]->size) 
            return;
        mpitracef("%d %d", i, bufsize);
        free(buffers[i]);
        buffers[i] = mpi_driver::newMpiBuffer(bufsize);
    }
};


mpi_driver::mpi_driver() 
{
    hostname = NULL;
};

mpi_driver::~mpi_driver()
{
    delete[] hostname;
    mpitracef("destruct %d", rank);
}

mpiBuffer * mpi_driver::newMpiBuffer()
{
    return newMpiBuffer(NOVO_MPI_MSG_LEN);
}

mpiBuffer * mpi_driver::newMpiBuffer(int size)
{
    mpiBuffer *b = (mpiBuffer*)malloc(sizeof(mpiBuffer) + size); 
    if(b == NULL)
    {
        fprintf(stderr, "%s:%s:%d  Memory Allocation Error\n", __FILE__, __func__, __LINE__);
        exit(-1);
    }
    b->size = size;
    return b;    
}


//Starts MPI Process and decides whether we are a master or slave
void mpi_driver::run_mpi(int argc, char** argv)
{
    int rc;
    int got_threading;
    rc = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &got_threading);
    mpitracef("mpi_driver::run_mpi() Threading %d got %d", MPI_THREAD_FUNNELED, got_threading);
    if (rc != MPI_SUCCESS) 
    {
        fprintf (stderr, "Error initialising %s, MPI::Init failed. Terminating.\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    mpitracel;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    mpitracel;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    if (num_tasks <= 1) 
    {
        fprintf (stderr, "Error initialising %s, no slave tasks. Terminating.\n"
                "You need to start at least two processes. eg. mpiexec -np 2 novoalignMPI ...\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    mpitracef("Rank %d, Num %d", rank, num_tasks);
    hostname = new char[MPI_MAX_PROCESSOR_NAME];
    (void)MPI_Get_processor_name(hostname, &l_hostname );
    if(rank != 0)
        run_slaves(argc, argv);
    else
        run_master(argc, argv);

    MPI_Finalize();

}
#define YIELD 1000000    // Initial yield is 1msec.
void mpi_driver::run_master(int argc, char** argv)
{
    static const int SENDBUFFERSPERSLAVE = 2;
    static const int RECVBUFFERS = num_tasks;
    int task, buf_idx;
    unsigned int tag;
    int block_num = 0;
    masterRequest = new msgbuffers;
    
    #define STDIOBUFSZ 64*1024
    char *stdioBuf = NULL;
    #ifndef _DEBUG
        setvbuf(stdout, stdioBuf = (char *)malloc(STDIOBUFSZ), _IOFBF, STDIOBUFSZ);
    #endif

    mpitracel;
    init_master(argc, argv);
    mpitracel;
    masterRequest->allocate(RECVBUFFERS, RECVBUFFERS + num_tasks*SENDBUFFERSPERSLAVE);
    eof_sent = new bool[num_tasks];
    wait_sync = new bool[num_tasks];
    for(int i = 0; i < num_tasks; i++)
        eof_sent[i] = false;
    mpitracel;
    int num_running;
#ifdef YIELD
    struct timespec sleep_for;
    sleep_for.tv_nsec = YIELD;
    sleep_for.tv_sec = 0;
#endif
    for(int i = 1; i < RECVBUFFERS; i++)
        MPI_Irecv(masterRequest->buffers[i]->message, masterRequest->buffers[i]->size, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &(masterRequest->requests)[i]);
    do 
    {
        num_running = num_tasks - 1;
        for(int i = 0; i < num_tasks; i++)
            wait_sync[i] = false;
        mpitracel;
        for(buf_idx = RECVBUFFERS; buf_idx < masterRequest->size ; buf_idx++) 
        {
            task = (buf_idx - RECVBUFFERS) % num_tasks;
            if(task == 0) 
                continue;
            send_message(task, buf_idx);
        }
        mpitracel;

        while(num_running > 0)
        {
            int flag = true, req_idx;
            MPI_Status status;
#ifndef YIELD
            // On some MPI's WaitAny sits in CPU Loops. 
            MPI_Waitany(masterRequest->size, masterRequest->requests, &buf_idx, &status);
#else
//            fprintf(stderr, "Testany %ld\n", sleep_for.tv_nsec);

            MPI_Testany(masterRequest->size, masterRequest->requests, &buf_idx, &flag, &status);
            if(!flag || buf_idx == MPI_UNDEFINED)
            {
                mpitracef("nanosleep R%d %ld", buf_idx % num_tasks, sleep_for.tv_nsec);
                nanosleep(&sleep_for, NULL);
                if(sleep_for.tv_nsec < 200000000)
                    sleep_for.tv_nsec += sleep_for.tv_nsec / 10;
                continue;
            }
            if(sleep_for.tv_nsec > 100000)
                sleep_for.tv_nsec -= sleep_for.tv_nsec / 2;
//            fprintf(stderr, "Success %ld\n", sleep_for.tv_nsec);

#endif
            mpitracef("testany R%d %d %d %d", buf_idx % num_tasks, buf_idx, num_tasks, num_running);
            if(buf_idx >= RECVBUFFERS)   // A Send completed
            {
                task = (buf_idx - RECVBUFFERS) % num_tasks;
                free_message(&(masterRequest->buffers[buf_idx]));
                send_message(task, buf_idx);
                continue;
            }
            // A receive completed
            mpitracef("Master Recvd %d from R%d", status.MPI_TAG, status.MPI_SOURCE);
            int i = buf_idx;
            int slave = status.MPI_SOURCE;
            if(status.MPI_TAG == TAG_MSGSIZE) 
            {
                unsigned int newsize = *(unsigned int *)(masterRequest->buffers[i]->message);
                mpitracef("From: R%d resize buf[%d] %d", slave, i, newsize);
                masterRequest->resize(i, newsize);
                MPI_Irecv(masterRequest->buffers[i]->message, masterRequest->buffers[i]->size, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &masterRequest->requests[i]);
                continue;
            }
            if(status.MPI_TAG == TAG_SYNCD) 
            {
                mpitracef("SYNCD from R%d Running %d Buf[%d]", slave, num_running, i);
                MPI_Irecv(masterRequest->buffers[i]->message, masterRequest->buffers[i]->size, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &masterRequest->requests[i]);
                num_running--;
                continue;
            }
            if(status.MPI_TAG == TAG_EOF) 
            {
                mpitracef("EOF from R%d Running %d Buf[%d]", slave, num_running, i);
                if(num_running > RECVBUFFERS)
                    MPI_Irecv(masterRequest->buffers[i]->message, masterRequest->buffers[i]->size, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &masterRequest->requests[i]);
                num_running--;
                continue;
            }
            (void)MPI_Get_count(&status, MPI_CHAR, &(masterRequest->buffers[i]->used));
            master_process_message(slave, masterRequest->buffers[buf_idx]->message, status.MPI_TAG, masterRequest->buffers[i]->used);
            MPI_Irecv(masterRequest->buffers[i]->message, masterRequest->buffers[i]->size, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &masterRequest->requests[i]);
        }
        if(!eof())
            synced();
    } while(!eof());
    mpitracef("Barrier");
    MPI_Barrier(MPI_COMM_WORLD);
    final_results();
    mpitracef("exit");
    if(stdioBuf != NULL) 
    {
        fflush(stdout);
        setvbuf(stdout, NULL, _IOFBF, 4*1024);
        free(stdioBuf);
    }

}

void mpi_driver::send_message(int task, int requestIdx)
{
    masterRequest->buffers[requestIdx] = next_message();
    if(masterRequest->buffers[requestIdx] != NULL)
    {
        mpitracef("Issend task R%d buf %d len %d", task, requestIdx, masterRequest->buffers[requestIdx]->used);
        MPI_Issend(masterRequest->buffers[requestIdx]->message, masterRequest->buffers[requestIdx]->used, MPI_CHAR, task, masterRequest->buffers[requestIdx]->tag, MPI_COMM_WORLD, &masterRequest->requests[requestIdx]);
        mpitracef("Issend task R%d buf %d len %d", task, requestIdx, masterRequest->buffers[requestIdx]->used);
    } 
    else if(eof()) 
    {
        if(!eof_sent[task]) 
        {
            mpitracef("send eof R%d", task);
            MPI_Issend((void *)eof_message, eof_len, dtConvert[eof_datatype], task, TAG_EOF, MPI_COMM_WORLD, &masterRequest->requests[requestIdx]);
            eof_sent[task] = true;
        }
    }
    else if(syncPoint())
    {
        if (!wait_sync[task]) 
        {
            mpitracef("send SYNC R%d", task);
            MPI_Issend((void *)eof_message, 0, dtConvert[eof_datatype], task, TAG_SYNC, MPI_COMM_WORLD, &masterRequest->requests[requestIdx]);
            wait_sync[task] = true;
        }
    }
}

void mpi_driver::init_master(int argc, char** argv) {};

bool mpi_driver::final_results() { return true; };

void mpi_driver::run_slaves(int argc, char** argv)
{
    mpitracel;
    init_slave(argc, argv);
    mpitracel;
    eof_received = false;
//    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Status status;
    while(true)
    {
        unsigned char eof_buf[eof_len];
        mpitracef("sR%d wait Probe()", rank);
        MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        mpitracef("sR%d probe success", rank);
        int tag = status.MPI_TAG;
        if(tag == TAG_EOF) 
        {
            mpitracef("sR%d eof recvd", rank);
            MPI_Recv(eof_buf, mpi_driver::eof_len,  dtConvert[mpi_driver::eof_datatype], 0, TAG_EOF, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            break;
        }
        mpitracef("sR%d",rank);
        int msg_len;
        (void)MPI_Get_count(&status, MPI_CHAR, &msg_len);
        process_message(tag, msg_len);
    }
    eof_slave();
    Barrier();
    final_slave();
    mpitracef("sR%d exit slave", rank);
}

void mpi_driver::init_slave(int argc, char** argv) {};

void mpi_driver::eof_slave() {};

void mpi_driver::final_slave() {};

int NV_MPI_Send(void* buf, int count, NV_Datatype datatype, int task, int tag)
{
    mpitracef("R%d %lx %d %d %x", mpidriver::rank, (long unsigned int)buf, count, datatype, dtConvert[datatype]);
    return MPI_Send(buf, count, dtConvert[datatype], task, tag, MPI_COMM_WORLD);
}
int NV_MPI_Recv(void * msg, int msg_len,  NV_Datatype datatype)
{
    mpitracef("R%d L:%d", mpi_driver::rank, msg_len);
    return MPI_Recv(msg, msg_len,  dtConvert[datatype], MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

int NV_MPI_Recv(void * msg, int msg_len,  NV_Datatype datatype, NV_MPI_STATUS status)
{
    mpitracef("R%d L:%d st:%lx", rank, msg_len, (long unsigned int)status);
    return MPI_Recv(msg, msg_len,  dtConvert[datatype], MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, (MPI_Status *)status);
}

int NV_MPI_MsgLen(int count, NV_Datatype datatype)
{
    MPI_Aint size;
    MPI_Aint lb;
    MPI_Type_get_extent(dtConvert[datatype], &lb, &size);
    mpitracef("R%d sz:%d n:%d", rank, (int)size, count);
    return size * count;
}

void Die()
{
    mpitracef("Abort[%d]", rank);
    MPI_Abort(MPI_COMM_WORLD, -1);
}

int NV_MPI_Bcast(void* buf, int count, NV_Datatype datatype)
{
    mpitracef("BCast[%d] %d ", rank, count);
    return MPI_Bcast(buf, count, dtConvert[datatype], 0, MPI_COMM_WORLD);
}

int Barrier()
{
    mpitracef("Barrier[%d]", rank);
    return MPI_Barrier(MPI_COMM_WORLD);
}

NV_MPI_STATUS NV_Status_Object() 
{
    return new MPI_Status; 
}

void deleteNV_Status_Object(NV_MPI_STATUS st)
{
    delete (MPI_Status *)st;
}

int NV_MPI_Get_count(NV_MPI_STATUS status, NV_Datatype datatype)
{
    int msg_len;
    MPI_Get_count((MPI_Status *)status, dtConvert[datatype], &msg_len);
    mpitracef("MPI_Get_count[%d] %d", rank, msg_len);
    return msg_len;
}

int mpi_driver::max_length() { return NOVO_MPI_MSG_LEN; };
