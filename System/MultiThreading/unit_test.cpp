#include <iostream>
#include <boost/thread/thread.hpp>
#include "producer_consumer.h"

void produce(ProducerConsumer<int> * Q);

void consume(ProducerConsumer<int> * Q);

int main()
{
    ProducerConsumer<int> * Q = new ProducerConsumer<int>(10);
    boost::thread * p_thr  = new boost::thread ( produce , Q );
    boost::thread * p_thr1 = new boost::thread ( produce , Q );
    boost::thread * p_thr2 = new boost::thread ( produce , Q );
    boost::thread * p_thr3 = new boost::thread ( produce , Q );
    boost::thread * c_thr  = new boost::thread ( consume , Q );
    boost::thread * c_thr1 = new boost::thread ( consume , Q );
    boost::thread * c_thr2 = new boost::thread ( consume , Q );
    boost::thread * c_thr3 = new boost::thread ( consume , Q );
    while(Q->size())
    {
        std::cout << Q->size() << std::endl;
        sleep(1);
    }
    p_thr  -> join ();
    p_thr1 -> join ();
    p_thr2 -> join ();
    p_thr3 -> join ();
    c_thr  -> join ();
    c_thr1 -> join ();
    c_thr2 -> join ();
    c_thr3 -> join ();
    std::cout << Q -> size () << std::endl;
    return 0;
}

void produce(ProducerConsumer<int> * Q)
{
    for ( int k(0)
        ; k < 10
        ; ++k
        )
    {
        Q -> produce ( k );
    }
}

void consume(ProducerConsumer<int> * Q)
{
    for ( int k(0)
        ; k < 10
        ; ++k
        )
    {
        int i = Q ->  consume (  );
        std::cout << "i=" << i << std::endl;
    }
}

