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
   bool sphere_line_contact(double a[4], double nodes[][3], int facets[][4] ); 
   bool contact_within_facet(double a[4], double nodes[][3], int facets[][4] );

   int facet_size;
  
   double Normal_facet_vector[3];
   double distance;
   bool debug_contact;
   double c[3];
   double b[3];
   double N_value;
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

  account1.debug_contact = false;

  bool contact;

/*  contact = account1.sphere_sphere_contact( sphere1, sphere2);

  if (contact)
    cout << "They are in contact " << endl;
  else if (!contact)
    cout << "They are not in contact " << endl;
*/

  contact = account1.sphere_facet_contact(sphere2, nodes, facets);

  if (contact)
  {
   contact =  account1.contact_within_facet(sphere2, nodes, facets);
   
   if (contact)
   {
     cout << "Contacted Facet " << endl;
   }

   else
   {
     contact = account1.sphere_line_contact(sphere2, nodes, facets);

     if (contact)
     {
       cout << "Contact line on facet " << endl;
     }

     else 
       cout << "Not in contact with facet " << endl;
   }
  }
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
  if (debug_contact)
    cout << endl << "Enter sphere_sphere_contact " << endl;

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
  if (debug_contact)
    cout << endl << "Enter Sphere Facet Contact" << endl;

  N_value = 0;
  double N_d;
  bool value;

// Distance from plane containing facet and sphere is less than the sphere radius

// A loop to check for sphere against all facets

  for (int i = 0; i < facet_size; i++ )
  {

// gets the two unit vectors that define the facet

    for (int j = 0; j < 3; j++)
    {
      c[j] = nodes[facets[i][0]][j] - nodes[facets[i][2]][j];
      b[j] = nodes[facets[i][1]][j] - nodes[facets[i][2]][j];

      if (debug_contact)
      {
        cout <<"Vector 1 in "<< j <<" direction = " << c[j] << endl;
        cout <<"Vector 2 in "<< j <<" direction = " << b[j] << endl;
      } 
    }

// determines the magnitude and direction of Normal vector to the plane made by the facet

    cross_product( c, b , Normal_facet_vector );

    if (debug_contact)
    {
      for (int j = 0; j < 3; j++)
        cout << " Normal facet vector in " << j << " direction = " << Normal_facet_vector[j] << endl;
    }

// converts the vector into a unit vector (direction)
  
    for (int j = 0; j < 3; j++)
      N_value += (Normal_facet_vector[j] * Normal_facet_vector[j]);
    N_value = sqrt( N_value );
    for (int j = 0; j < 3; j++)
      Normal_facet_vector[j] = Normal_facet_vector[j] / N_value;

    if (debug_contact)
    {
      for (int j= 0; j < 3; j++)
        cout << "Unit Normal facet vector in " << j << " direction = " << Normal_facet_vector[j] << endl;
    }

// determines the value of Nd by dotting Normal vector with a point on the facet 

    N_d = Normal_facet_vector[0] * nodes[facets[i][0]][0] + Normal_facet_vector[1] * nodes[facets[i][0]][1] + Normal_facet_vector[2] * nodes[facets[i][0]][2];

    if (debug_contact)
      cout << "Nd value for the facet plane = " << N_d << endl;

// determines the distance from the facet plane to the sphere center location

    distance = dot_product(Normal_facet_vector, a) + N_d;

    if (debug_contact)
      cout << "Distance from the sphere to the facet plane = " << distance << endl;

// if the distance from sphere center to plane is less than radius, check if the contact point lies within the boundaries of the facet

   if (sqrt(distance*distance < a[3] ) )
   {
     value = true;  
   }
    else
      value = false;
  }
  
return value;
     
}


bool math::sphere_line_contact(double a[4], double nodes[][3], int facets[][4] )
{
  if (debug_contact)
    cout << endl << "Enter Sphere Line Contact " << endl;

  double vector_ca[3], vector_cb[3], length_l, length_L = 0, abc[3];
  distance = 0;
  double Area_abc; 
  bool value;

// loop to check contact on every facet

  for (int i = 0; i < facet_size; i++ )
  {
    if (debug_contact)
      cout << "Entering Facet Number #" << i << endl;

// loop to check all three lines on each facet

    for (int j = 0; j < 3; j++)
    {

// loop to obtain vector between two nodes and vector from node to sphere center

      for (int k = 0; k < 3; k++)
      {
        if (j == 0)
        {
          vector_ca[k] = nodes[facets[i][0]][k] - nodes[facets[i][1]][k];
          vector_cb[k] = a[k] - nodes[facets[i][0]][k];
  
        }
        else if (j == 1)
        {
          vector_ca[k] = nodes[facets[i][1]][k] - nodes[facets[i][2]][k];
          vector_cb[k] = a[k] - nodes[facets[i][1]][k];
        }
        else 
        {
          vector_ca[k] = nodes[facets[i][0]][k] - nodes[facets[i][2]][k];
          vector_cb[k] = a[k] - nodes[facets[i][0]][k];
        }

        if (debug_contact)
        {
          cout << "Vector CA in the " << k << " direction = " << vector_ca[k] << endl;
          cout << "Vector CB in the " << k << " direction = " << vector_cb[k] << endl;
        }
      }

// the length between the node and the intersection of the perpendicular distance to sphere center

      length_l = dot_product(vector_ca, vector_cb);
 
      if (debug_contact)
        cout << "Length l = " << length_l << endl;

// the length of the line (between two nodes)

      length_L = sqrt( vector_ca[0] * vector_ca[0] + vector_ca[1] * vector_ca[1] + vector_ca[2] * vector_ca[2] );

      if (debug_contact)
        cout << "Length L = " << length_L << endl;

// to determine the distance between line and the sphere center : d = 2 * triangle made by sphere center / length of line

      cross_product( vector_ca, vector_cb, abc);

      Area_abc = sqrt( abc[0] * abc[0] + abc[1] * abc[1] + abc[2] * abc[2] );

      distance = 2 * Area_abc / length_L; 

      if (debug_contact)
        cout << "Distance = " << distance << endl << endl;
  
// if the point of the perpendicular line to sphere center is between the two nodes 

      if ( (abs(length_l) >= 0 ) && (abs(length_l) <= length_L) )
      {

// if the perpendicular distance to line is less than radius of sphere, then contact

        if (distance < a[i] )
        {
          value = true;
          break;
        }

        else
          value = false;
      }

      else
        value = false;

    }

  }

return value;

}
       
bool math::contact_within_facet(double a[4], double nodes[][3], int facets[][4] )
{
  if (debug_contact)
    cout << endl << "Enter Contact Within Facet " << endl;

  bool value;

  double A[3], abc[3], B[3], areaA = 0, areaB = 0, areaABC = 0, alpha = 0, beta = 0, gamma = 0;

  double point_inside[3];

// translates the sphere center to a point on the plane

  for (int k = 0; k < 3; k++)
    point_inside[k] = a[k] - distance * Normal_facet_vector[k];

  if (debug_contact)
  {
    for (int j = 0; j < 3; j++)
      cout << "Location of point on plane for " << j << " direction = " << point_inside[j] << endl;
  }

// determines the area for the sphere point on plane with a facet node 

   cross_product( c , point_inside, A);

   for (int k = 0; k < 3; k++ )
     areaA += (A[k] * A[k] );
 
   areaA = sqrt( areaA );

   cross_product( b , point_inside, B);

   for (int k = 0; k < 3; k++ )
     areaB += (B[k] * B[k] );
 
   areaB = sqrt( areaB );
 
   for (int k = 0; k < 3; k++)
     areaABC += (Normal_facet_vector[k] * N_value * Normal_facet_vector[k] * N_value );

   areaABC = sqrt( areaABC );

// alpha , beta, gamma based upon Barycentric Coordinates

   alpha = areaA / areaABC ; 
   beta = areaB / areaABC ;
   gamma = 1 - alpha - beta;

// if point lies within the triangle, all three coordinates are less than or equal to one
// one of the coordinates equals 0 if the point lies on the triangle boundary
// one of the coordinates equal 1 if the point lies on a vertex

   if ( (alpha <= 1 ) && (beta <= 1 ) && (gamma >= 0 ) )
     value = true;
   else
   {
     value = false;
   }

  return value;
} 

