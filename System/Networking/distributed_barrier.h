#ifndef DISTRIBUTED_BARRIER_H
#define DISTRIBUTED_BARRIER_H

#include "networking.h"
#include "producer_consumer.h"

void distributed_barrier_server_worker ( int k 
                                       , int n
                                       , Server * s 
                                       , ProducerConsumer < int > * Q
                                       )
{
    Message m;
    s . Receive ( m , k );
    Q . produce ( k );
    if ( Q . size () == n )
    {
        // broadcast unlock signal to clients
        for ( int j(0)
            ; j < n
            ; ++j
            )
        {
            s . Send ( m , j );
        }
    }
}

class DistributedBarrier_server
{
    ProducerConsumer < int > * Q;
    Server * s;
    DistributedBarrier_server ( long port
                              , int num_clients
                              )
    {
        Q = new ProducerConsumer < int > ( num_clients );
        s = new Server;
        s -> Create ( port , num_clients );
        for ( int k(0)
            ; k < num_clients
            ; ++k
            )
        {
            boost::thread thr = new boost::thread ( distributed_barrier_server_worker 
                                                  , k
                                                  , num_clients
                                                  , s
                                                  , Q
                                                  );
        }
    }
};

class DistributedBarrier_client
{
    Client c;
    DistributedBarrier_client ( std::string host
                              , int port
                              )
    {
        c . Connect ( host , port );
    }
    void Block ()
    {
        Message m;
        c . Receive (m);
    }
    void Signal ()
    {
        Message m;
        Test t;
        m << t;
        c . Send (m);
    }
};


#endif

