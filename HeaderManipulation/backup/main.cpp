#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <stack>
#include <deque>

using namespace std;

void Tokenize(string const & str, vector<string>& tokens, const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    string::size_type pos0 = str.find_first_of(delimiters, 0);
    while ( string::npos != pos0 )
    {
        if ( pos0 > str.size() )
        {

        }
        else
        {
            if ( pos0 < lastPos )
            {
                tokens.push_back(str.substr(pos0,1));
            }
        }
        pos0 = str.find_first_of(delimiters,pos0+1);
    }
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        string::size_type p_lastPos ( lastPos ) , p_pos;
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        while(1)
        {
            if ( pos >= str.size() )
            {

            }
            else
            {
                tokens.push_back(str.substr(pos,1));
            }
            p_pos = str.find_first_of(delimiters, pos+1);
            if ( string::npos == pos 
              || string::npos == p_pos
              || string::npos == p_lastPos
              || p_pos+1 > lastPos 
               )
            {
                break;
            }
            pos = p_pos;
        }
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

bool operator == ( std::string const & s , char const * ch )
{
    int count = 0;
    char const * c ( &ch[0] );
    while ( *c != '\0' )
    {
        if ( count >= s.size() ) return false;
        if ( s[count] != *c ) return false;
        ++count;
        ++c;
    }
    if ( count != s.size() ) return false;
    return true;
}

bool operator == ( std::string & s , char const * ch )
{
    int count = 0;
    char const * c ( &ch[0] );
    while ( *c != '\0' )
    {
        if ( count >= s.size() ) return false;
        if ( s[count] != *c ) return false;
        ++count;
        ++c;
    }
    if ( count != s.size() ) return false;
    return true;
}

bool operator == ( std::string const & s , char * ch )
{
    int count = 0;
    char const * c ( &ch[0] );
    while ( *c != '\0' )
    {
        if ( count >= s.size() ) return false;
        if ( s[count] != *c ) return false;
        ++count;
        ++c;
    }
    if ( count != s.size() ) return false;
    return true;
}

bool operator == ( std::string & s , char * ch )
{
    int count = 0;
    char const * c ( &ch[0] );
    while ( *c != '\0' )
    {
        if ( count >= s.size() ) return false;
        if ( s[count] != *c ) return false;
        ++count;
        ++c;
    }
    if ( count != s.size() ) return false;
    return true;
}

enum type { OPERAND           = 1
          , OPERATOR          = 2
          , FUNCTION          = 3
          , OPEN_PARENTHESIS  = 4
          , CLOSE_PARENTHESIS = 5
          , COMMAS            = 6
          , SPACE             = 100
          };

struct token
{
    string id;
    type t;
    token ( string const & _id )
    : id ( _id )
    , t ( get_type() )
    {

    }
    bool is_function() const
    {
        if ( t == FUNCTION ) return true;
        return false;
    }
    type get_type() const
    {
        if ( id == "," ) return COMMAS;
        if ( id == " " 
           , id == "\t"
           , id == "\n"
           ) 
        {
            return SPACE;
        }
        if ( id == "(" ) return OPEN_PARENTHESIS;
        if ( id == ")" ) return CLOSE_PARENTHESIS;
        if ( id == "cos" 
          || id == "acos"
          || id == "sin"
          || id == "asin"
          || id == "tan"
          || id == "atan"
          || id == "sqrt"
          || id == "exp"
          || id == "pow"
          || id == "abs"
          || id == "negate"
          || id == "int"
          || id == "if"
           )
        {
            return FUNCTION;
        }
        if ( id == "+"
          || id == "-"
          || id == "*"
          || id == "/"
          || id == "^"
          || id == "%"
          || id == "="
           )
        {
            return OPERATOR;
        }
        return OPERAND;
    }
    int get_presedence() const
    {
        if ( t == "=" )
        {
            return -1;
        }
        if ( t == COMMAS  )
        {
            return 1;
        }
        if ( id == "%" )
        {
            return 2;
        }
        if ( id == "^" )
        {
            return 4;
        }
        if ( id == "+" || id == "-" )
        {
            return 5;
        }
        if ( id == "*" || id == "/" )
        {
            return 10;
        }
        if ( id == "(" || id == ")" || t == FUNCTION )
        {
            return 0;
        }
        return -1000;
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// While there are tokens to be read:
// 
//     Read a token.
//     If the token is a number, then add it to the output queue.
//     If the token is a function token, then push it onto the stack.
//     If the token is a function argument separator (e.g., a comma):
// 
//         Until the token at the top of the stack is a left parenthesis, pop operators off the stack onto the output queue. 
//         If no left parentheses are encountered, either the separator was misplaced or parentheses were mismatched.
// 
//     If the token is an operator, o1, then:
// 
//         while there is an operator token, o2, at the top of the stack, and
// 
//                 either o1 is left-associative and its precedence is less than or equal to that of o2,
//                 or o1 is right-associative and its precedence is less than that of o2,
// 
//             pop o2 off the stack, onto the output queue;
// 
//         push o1 onto the stack.
// 
//     If the token is a left parenthesis, then push it onto the stack.
//     If the token is a right parenthesis:
// 
//         Until the token at the top of the stack is a left parenthesis, pop operators off the stack onto the output queue.
//         Pop the left parenthesis from the stack, but not onto the output queue.
//         If the token at the top of the stack is a function token, pop it onto the output queue.
//         If the stack runs out without finding a left parenthesis, then there are mismatched parentheses.
// 
// When there are no more tokens to read:
// 
//     While there are still operator tokens in the stack:
// 
//         If the operator token on the top of the stack is a parenthesis, then there are mismatched parentheses.
//         Pop the operator onto the output queue.
// 
// Exit.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void parse ( string const & str )
{

    std::cout << str << std::endl;

    vector < string > t;

    Tokenize ( str , t , " +-*/%()," );

    vector < token > v;

    for ( std::size_t i(0)
        ; i < t.size()
        ; ++i
        )
    {
        v . push_back ( token ( t[i] ) );
    }

    std::stack < token > S;

    std::deque < token > Q;

    for ( std::size_t i(0)
        ; i < v.size()
        ; ++i
        )
    {
        if ( v[i].get_type() == FUNCTION )
        {
            S . push ( v[i] );
        }
        else
        if ( v[i].get_type() == OPEN_PARENTHESIS )
        {
            S . push ( v[i] );
        }
        else
        if ( v[i].get_type() == CLOSE_PARENTHESIS )
        {
            bool found_matching_parenthesis = false;
            while ( !S.empty() )
            {
                token tmp = S . top ();
                if ( tmp.get_type() == OPEN_PARENTHESIS )
                {
                    S . pop ();
                    found_matching_parenthesis = true;
                    break;
                }
                else
                {
                    S . pop ();
                    Q . push_front ( tmp );
                }
            }
            if ( found_matching_parenthesis == false )
            {
                std::cout << "missing parenthesis." << std::endl;
                exit(1);
            }
            if ( !S.empty() )
            {
                token fun = S . top ();
                if ( fun . is_function() )
                {
                    S . pop ();
                    Q . push_front ( fun );
                }
            }
        }
        else
        if ( v[i].get_type() == OPERATOR || v[i].get_type() == COMMAS )
        {
            while ( !S.empty() )
            {
                token tmp = S . top ();
                if ( v[i].get_presedence() > tmp.get_presedence() )
                {
                    break;
                }
                S . pop ();
                Q . push_front ( tmp );
            }
            S . push ( v[i] );
        }
        else
        if ( v[i].get_type() == OPERAND )
        {
            Q . push_front ( v[i] );
        }
    }

    while ( !S.empty() )
    {
        token tmp = S . top ();
        S . pop ();
        Q . push_front ( tmp );
    }

    for ( std::deque < token > :: iterator it ( Q . begin () )
        ; it != Q . end ()
        ; ++it
        )
    {
        std::cout << it->id << std::endl;
    }

}

int main(int argc, char *argv[])
{
    string s = "neg(3/(cos(5+3*9)+neg(1*PATH7)))";

    parse ( s );

    return 0;
}


