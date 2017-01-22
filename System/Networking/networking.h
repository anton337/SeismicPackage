#ifndef NETWORKING_H
#define NETWORKING_H

/* A simple server in the internet domain using TCP
   The port number is passed as an argument */

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

const long chunk_size = 512;//4*5000;

struct Test
{
    float x;
    float y;
    float z;
    Test (){}
    Test ( float _x
         , float _y 
         , float _z
         )
    : x ( _x )
    , y ( _y )
    , z ( _z )
    {

    }
};

struct Chunk
{
    char buffer[chunk_size];
    Chunk()
    {
        bzero(buffer,chunk_size);
    }
};

struct Message
{
    std::vector < Chunk > data;
};

Message & operator << ( Message & m , Test const & t )
{
    Chunk c;
    float * x ( reinterpret_cast < float * > ( & c . buffer [ 0 ] ) );
    *x = t . x;
    float * y ( reinterpret_cast < float * > ( & c . buffer [ 4 ] ) );
    *y = t . y;
    float * z ( reinterpret_cast < float * > ( & c . buffer [ 8 ] ) );
    *z = t . z;
    m . data . push_back ( c );
};

Message & operator >> ( Message & m , Test & t )
{
    Chunk c = m . data[0];
    float * x ( reinterpret_cast < float * > ( & c . buffer [ 0 ] ) );
    t . x = *x;
    float * y ( reinterpret_cast < float * > ( & c . buffer [ 4 ] ) );
    t . y = *y;
    float * z ( reinterpret_cast < float * > ( & c . buffer [ 8 ] ) );
    t . z = *z;
    return m;
};

void error(const char *msg)
{
    perror(msg);
    exit(1);
}

struct Server
{
    int sockfd, portno;
    std::vector < int > newsockfd;
    void Create ( int port , int num_client )
    {
         std::cout << "CREATE " << num_client << std::endl;
         struct sockaddr_in serv_addr;
         sockfd = socket(AF_INET, SOCK_STREAM, 0);
         if (sockfd < 0) 
            error("ERROR opening socket");
         bzero((char *) &serv_addr, sizeof(serv_addr));
         portno = port;
         serv_addr.sin_family = AF_INET;
         serv_addr.sin_addr.s_addr = INADDR_ANY;
         serv_addr.sin_port = htons(portno);
         int yes = 1;
         if ( setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) == -1 )
         {
             error("setsockopt");
         }
         if (bind(sockfd, (struct sockaddr *) &serv_addr,
                  sizeof(serv_addr)) < 0) 
                  error("ERROR on binding");
         listen(sockfd,50);
         for ( int k(0)
             ; k < num_client
             ; ++k
             )
         {
            std::cout << "ACCEPT " << k << std::endl;
            socklen_t clilen;
            struct sockaddr_in cli_addr;
            clilen = sizeof(cli_addr);
            int newfd = accept(sockfd, 
                              (struct sockaddr *) &cli_addr, 
                              &clilen
                              ); 
            std::cout << "ACCEPTed " << k << std::endl;
            if (newfd < 0) 
                 error("ERROR on accept");
            newsockfd . push_back ( newfd );
            std::cout << "newsockfd size : " << newsockfd.size() << std::endl;
         }
    }
    void Send ( Message const & m , int k , int num_bytes = chunk_size )
    {
         if ( k < newsockfd.size () )
         {
         int num = num_bytes / chunk_size;
         if ( num != m . data . size () ) {std::cout << num << " " << m . data . size () << std::endl; error ( "server message size mismatch." );}
         char buffer[chunk_size];
         int n;
         for ( int i(0)
             ; i < num//m . data . size ()
             ; ++i
             )
         {
            n = write(newsockfd[k],m . data[i] . buffer,chunk_size);
            if (n < 0) error("ERROR writing to socket");
         }
         }
         else
         {
            //error("attempting to send before socket has been accepted.");
            sleep(1);
         }
    }
    void Receive ( Message & m , int k , int num_bytes = chunk_size )
    {
         if ( k < newsockfd.size () )
         {
         int num = num_bytes / chunk_size;
         char buffer[chunk_size];
         int n;
         //while ( 1 )
         for ( int i(0)
             ; i < num//m . data . size ()
             ; ++i
             )
         {
            bzero(&buffer[0],chunk_size);
            n = read(newsockfd[k],&buffer[0],chunk_size);
            //sleep(1);
            if (n < 0) error("ERROR reading from socket");
            Chunk c;
            memcpy ( c . buffer , buffer , chunk_size );
            m . data . push_back (c);
            usleep(5000);
         }
         }
         else
         {
            //error("attempting to receive before socket has been accepted.");
            sleep(1);
         }
    }
    void Destroy ()
    {
         for ( int k(0)
             ; k < newsockfd.size()
             ; ++k
             )
         close(newsockfd[k]);
         close(sockfd);
    }
};

struct Client
{
    int sockfd, portno;
    char buffer[chunk_size];
    void Connect(std::string host,int port)
    {
        struct sockaddr_in serv_addr;
        struct hostent *server;
        portno = port;
        sockfd = socket(AF_INET, SOCK_STREAM, 0);
        if (sockfd < 0) 
            error("ERROR opening socket");
        server = gethostbyname(host.c_str());
        if (server == NULL) {
            fprintf(stderr,"ERROR, no such host %s\n",host.c_str());
            exit(0);
        }
        bzero((char *) &serv_addr, sizeof(serv_addr));
        serv_addr.sin_family = AF_INET;
        bcopy((char *)server->h_addr, 
             (char *)&serv_addr.sin_addr.s_addr,
             server->h_length);
        serv_addr.sin_port = htons(portno);
        if (connect(sockfd,(struct sockaddr *) &serv_addr,sizeof(serv_addr)) < 0) 
            error("ERROR connecting");
    }
    void Send ( Message const & m , int num_bytes = chunk_size )
    {
         int num = num_bytes / chunk_size;
         if ( num != m . data . size () ) {std::cout << num << " " << m . data . size () << std::endl; error ( "client message size mismatch." );}
         int n;
         for ( int i(0)
             ; i < num /*m . data . size ()*/
             ; ++i
             )
         {
            n = write(sockfd,m . data[i] . buffer,chunk_size);
            if (n < 0) error("ERROR writing to socket");
         }
    }
    void Receive ( Message & m , int num_bytes = chunk_size )
    {
         int num = num_bytes / chunk_size;
         int n;
         //while ( 1 )
         for ( int i(0)
             ; i < num /*m . data . size ()*/
             ; ++i
             )
         {
            bzero(buffer,chunk_size);
            n = read(sockfd,buffer,chunk_size);
            if (n < 0) error("ERROR reading from socket");
            Chunk c;
            memcpy ( c . buffer , buffer , chunk_size );
            m . data . push_back (c);
         }
    }
    void Destroy ()
    {
        close(sockfd);
    }
};


#endif

