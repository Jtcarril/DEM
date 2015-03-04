#include<iostream>
#include<cmath>
#include<cstdlib>
#include<cstddef>

using namespace std;

class math
{

  public:
   double dot_product(double a[], double b[]);
   void cross_product(double a[3], double b[3], double c[3]);
   bool sphere_sphere_contact(double a[4], double b[4]);
   bool sphere_facet_contact(double a[4], double nodes[][3], int facets[][4] );
   bool sphere_line_contact(double a[4], double nodes_1[3], double nodes_2[3], double nodes_3[3]  );

   int facet_size;
};

int main()
{
  double sphere1[4] = {0,0,1,1};
  double sphere2[4] = {0.6,0.6,0,1};
  math account1;
  account1.facet_size = 1;
  double nodes[3][3];
  int facets[account1.facet_size][4];
  
  nodes[0][0] = 0;
  nodes[0][1] = 0;
  nodes[0][2] = 0;

  nodes[1][0] = 1;
  nodes[1][1] = 0;
  nodes[1][2] = 0;
 
  nodes[2][0] = 0;
  nodes[2][1] = 1;
  nodes[2][2] = 0;

  facets[0][0] = 0;
  facets[0][1] = 1;
  facets[0][2] = 2;
  facets[0][3] = 1;

  bool contact;

  contact = account1.sphere_sphere_contact( sphere1, sphere2);

  if (contact)
    cout << "They are in contact " << endl;
  else if (!contact)
    cout << "They are not in contact " << endl;

  contact = account1.sphere_facet_contact(sphere2, nodes, facets);

  if (contact)
    cout << "In contact with facet " << endl;
  else if (!contact)
    cout << "Not in contact with facet " << endl;

  return 0;

}

double math::dot_product(double a[], double b[])
{
  double value;
  
  value = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

  return(value);
}
 
void math::cross_product(double a[3], double b[3], double c[3])
{
  c[0] = a[1] * b[2] - b[1] * a[2];
  c[1] = -1 * ( a[0] * b[2] - b[0] * a[2] );
  c[2] = a[0] * b[1] - b[0] * a[1];
}  

bool math::sphere_sphere_contact(double a[4], double b[4])
{
  bool value;
  double distance = 0;

  for (int i = 0; i < 3; i++)
    distance += (b[i] - a[i]) * (b[i] - a[i]); 
  distance = sqrt( distance);

  if (distance < (a[3] + b[3] ) )
    value = true;
  else if (distance >= (a[3] + b[3] ) )
    value = false;
  else
  {
    cout << "Unknown Error[0]" << endl;
    return 0;
  }

  return value;
}

bool math::sphere_facet_contact(double a[4], double nodes[][3], int facets[][4])
{
  double N[3];
  double N_value = 0;
  double N_d;
  double distance;
  bool value;
  double c[3];
  double b[3];

  double point_inside[3];
  double A[3], abc[3], B[3], areaA = 0, areaB = 0, areaABC = 0, alpha = 0, beta = 0, gamma = 0;

  for (int i = 0; i < facet_size; i++ )
  {
    for (int j = 0; j < 3; j++)
    {
      c[j] = nodes[facets[i][0]][j] - nodes[facets[i][2]][j];
      b[j] = nodes[facets[i][1]][j] - nodes[facets[i][2]][j];

    }

    cross_product( c, b , N );
  
    for (int j = 0; j < 3; j++)
      N_value += (N[j] * N[j]);
    N_value = sqrt( N_value );
    for (int j = 0; j < 3; j++)
      N[j] = N[j] / N_value;
    N_d = N[0] * nodes[facets[i][0]][0] + N[1] * nodes[facets[i][0]][1] + N[2] * nodes[facets[i][0]][2];
    distance = dot_product(N, a) + N_d;

    if (sqrt(distance*distance) < a[3] )
    {
       for (int k = 0; k < 3; k++)
         point_inside[i] = a[i] - distance * N[i];

       cross_product( c , point_inside, A);

       for (int k = 0; k < 3; k++ )
         areaA += (A[k] * A[k] );
 
       areaA = sqrt( areaA );

       cross_product( b , point_inside, B);

       for (int k = 0; k < 3; k++ )
         areaB += (B[k] * B[k] );
 
       areaB = sqrt( areaB );
 
       for (int k = 0; k < 3; k++)
         areaABC += (N[k] * N_value * N[k] * N_value );

       areaABC = sqrt( areaABC );

       alpha = areaA / areaABC ; 
       beta = areaB / areaABC ;
       gamma = 1 - alpha - beta;

       if ( (alpha <= 1 ) && (beta <= 1 ) && (gamma >= 0 ) )
         value = true;
       else
       {
         double temp_1[3], temp_2[3], temp_3[3];
         for (int m = 0; m < 3; m++)
         {
           temp_1[m] = nodes[facets[i][0]][m];
           temp_2[m] = nodes[facets[i][1]][m];
           temp_3[m] = nodes[facets[i][2]][m];
         }
         cout << "Calling sphere_line " << endl;
         value =  sphere_line_contact(a, temp_1, temp_2, temp_3); 
       }
    }
    else
      value = false;
  }
  
return value;
     
}
bool math::sphere_line_contact(double a[4], double nodes_1[3], double nodes_2[3], double nodes_3[3]  )
{
  double vector_ca[3], vector_cb[3], vector_L[3], length_L = 0;
  double distance = 0;
  double  temp_nodes[3][3];

  for (int i = 0; i < 3; i++)
  {
    temp_nodes[0][i] = nodes_1[i];
    temp_nodes[1][i] = nodes_2[i];
    temp_nodes[2][i] = nodes_3[i];
  }

  for (int j = 0; j < 2; j++ )
  {

  for (int i = 0; i < 3; i++)
  {
    vector_ca[i] = a[i] - temp_nodes[j][i];
    vector_cb[i] = temp_nodes[j+1][i] - temp_nodes[j][i];
    length_L += (vector_ca[i] * vector_ca[i] );
  }

  length_L = sqrt( length_L );

  cross_product ( vector_ca, vector_cb, vector_L );

  for (int i = 0; i < 3; i++)
  {
    distance += (vector_L[i] * vector_L[i] );
  }

  distance = sqrt( distance);

  distance = 2 * distance / length_L ; 

  cout << "Distance = " << distance << endl;
  double length_l = dot_product( vector_ca, vector_cb );

  if ( (length_l >= 0 ) && (length_l <= length_L ) )
  {
    if (a[3] > distance )
      return true;
  }
  }

  for (int i = 0; i < 3; i++)
  {
    vector_ca[i] = a[i] - temp_nodes[3][i];
    vector_cb[i] = temp_nodes[0][i] - temp_nodes[2][i];
    length_L += (vector_ca[i] * vector_ca[i] );
  }

  length_L = sqrt( length_L );

  cross_product ( vector_ca, vector_cb, vector_L );

  for (int i = 0; i < 3; i++)
  {
    distance += (vector_L[i] * vector_L[i] );
  }

  distance = sqrt( distance);

  distance = 2 * distance / length_L ; 

  cout << "Distance = " << distance << endl;

  double length_l = dot_product( vector_ca, vector_cb );

  if ( (length_l >= 0 ) && (length_l <= length_L  ) )
  {
    if (a[3] > distance )
      return true;
    else 
      return false; 
  }
  else 
    return false;
}
  
