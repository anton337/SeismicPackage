#include "networking.h"
#include <string>
#include <boost/thread/thread.hpp>

void server_func (int port)
{
    Server s;
    s . Create(port,2);
    for ( int k(0)
        ; k < 2
        ; ++k
        )
    {
        Message m0;
        s . Receive(m0,k);
        Test o;
        m0 >> o;
        std::cout << o.x << " " << o.y << " " << o.z << std::endl;
        Message m1;
        Test t(5,6,7);
        m1 << t;
        s . Send(m1,k);
    }
    s . Destroy();
}

void client_func (int port)
{
    Client c;
    c . Connect("localhost",port);
    Message m0;
    Test t(1,2,3);
    m0 << t;
    m0 << t;
    c . Send(m0);
    Message m1;
    c . Receive(m1);
    Test o;
    m1 >> o;
    std::cout << o.x << " " << o.y << " " << o.z << std::endl;
    c . Destroy();
}

int main(int argc , char * argv[])
{
    if (argc != 2)
    {
        std::cout << "please specify port" << std::endl;
        exit(1);
    }
    int port = atoi(argv[1]);
    boost::thread * thr1 = new boost::thread ( server_func , port );
    boost::thread * thr2 = new boost::thread ( client_func , port );
    boost::thread * thr3 = new boost::thread ( client_func , port );
    thr1 -> join ();
    thr2 -> join ();
    thr3 -> join ();
    return 0;
}

