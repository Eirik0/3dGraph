package main;

public class Graph3d {
	const float pi = 3.14159265;
	const float two_pi = 2 * pi;
	//array used by main and is_it_closer function
	float distance_to_the_plane[1921][1201]; 
	//what_system   1 = spherical 2 = rectangular 4 = polar
	//where the graph and interface are drawn
	BITMAP *buffer;

	//for interface
	BITMAP *bar; 
	BITMAP *box; 
	BITMAP *rdra;
	BITMAP *rdra2;
	BITMAP *rota;   
	BITMAP *rota2;  

	//the equation for spherical coordinates
	float what_rho_is(float ph, float th, int which_type, float how_large)
	{
	  
	   float rho, theta, phi, returned_var;
	   
	   if(which_type == 1)
	   {
	      theta = th;
	      phi = ph;
	      /*-(rho)-*/
	      
	         rho = 1-sqrt(tan(sin(phi)));
	         
	      /*-(rho)-*/
	      returned_var=rho;
	   }
	   
	   if (which_type == 2)
	   {
	      phi = ph;
	      rho = (th / two_pi) * (750/(how_large+1));
	      /*-(theta)-*/
	      
	         theta = phi-pi*sin (rho);
	      
	      /*-(theta)-*/
	      returned_var=theta;
	   }
	   
	   if(which_type == 3)
	   {
	      theta = th;
	      rho = (ph / pi) * (750/(how_large+1));
	      /*-(phi)-*/
	      
	         phi = 1-rho*sin(theta);
	      
	      /*-(phi)-*/
	      returned_var=phi;
	   }
	  
	  
	   return returned_var;
	}



	//the equation for rectangular coordinates
	float what_z_is(float x1, float y1, int which_type)
	{
	   
	   float x, y, z, returned_var;
	   
	   if(which_type == 1)
	   {
	     x = x1;
	     y = y1;
	     /*-(z)-*/  
	      
	        z = sin(sqrt(x*y));
	        
	     /*-(z)-*/ 
	     returned_var = z;
	   } 
	   
	   if(which_type == 2)
	   {
	     x = x1;
	     z = y1;
	     /*-(y)-*/  
	      
	        y = x*x*x -  z*z;
	        
	     /*-(y)-*/ 
	     returned_var = y;
	   } 
	   
	   if(which_type == 3)
	   {
	     z = x1;
	     y = y1;
	     /*-(x)-*/  
	      
	        x = z + y*y -sin(z);
	        
	     /*-(x)-*/ 
	     returned_var = x;
	   }   
	   
	   return returned_var;
	}



	//the equation in polar
	float what_polar_z_is(float phi_va, float theta_va, int which_type,float how_large)
	{
	   float z, r, theta, returned_var;
	   
	   if(which_type == 1)
	   {
	      theta = theta_va;
	      r = phi_va;
	      /*-(z)-*/ 
	      
	         z = sqrt(r);
	         
	      /*-(z)-*/ 
	      returned_var = z;
	   }
	   
	   if(which_type == 2)
	   {
	      theta = theta_va;
	      z = phi_va; //(phi_va / pi) * (750/(how_large)) -((750/(how_large))/2);
	      /*-(r)-*/ 
	      
	         r = sin (z - theta);
	         
	      /*-(r)-*/ 
	      returned_var = r;
	   }
	   
	    if(which_type == 3)
	   {
	      r = phi_va;//(phi_va/ pi ) * (1500/how_large);
	      z = theta_va; //(theta_va / two_pi) * (750/(how_large)) -((750/(how_large))/2);
	      /*-(r)-*/ 
	      
	         theta = z-sin(cos(1-r));
	         
	      /*-(r)-*/ 
	      returned_var = theta;
	   }
	   
	   return returned_var;
	}




	//parametric equations
	void parametric_xyz(float& x_uv, float& y_uv, float& z_uv, float phi_u, float theta_v)
	{     
	   float u_min = 0;
	   float u_max = 2*pi;
	   float v_min = 0;
	   float v_max = 2*pi;
	   
	   float u, v;
	   
	   //change parameters
	   u = (phi_u / pi) * (u_max - u_min) + u_min;
	   v = (theta_v / two_pi) * (v_max - v_min) + v_min;
	   
	   /*-(parametric)-*/ 
	   
	    x_uv =  cos(v);
	   
	    y_uv =  sin(v)*cos(u);
	   
	    z_uv =  sin(u);

	   
	   
	   /*-(parametric)-*/ 
	   
	}



	void e_color(int &r, int &g, int &b, float col)
	{
	  if (col >=0 && col < 0.1)
	  {
	     r = 255;
	     g = int(145.0 * (col / 0.1));
	     b = int(30.0 * (col / 0.1));
	  }
	  else if (col >=0.1 && col < 0.2)
	  {
	     r = 255;
	     g = int(145.0 + 110.0 * ((col - 0.1) / 0.1));
	     b = int(30.0 + 20.0 * ((col - 0.1) / 0.1));
	  }
	  else if (col >=0.2 && col < 0.3)
	  {
	     r = int(255.0 - 80.0 * ((col - 0.2) / 0.1) );
	     g = 255;
	     b = int(50.0 - 50.0 * ((col - 0.2) / 0.1));
	  }
	  else if (col >=0.3 && col < 0.4)
	  {
	     r = int(175.0 - 175.0 * ((col - 0.3) / 0.1) );
	     g = 255;
	     b = 0;
	  }
	  else if (col >=0.4 && col < 0.5)
	  {
	     r = 0;
	     g = int(255.0 - 80.0 * ((col - 0.4) / 0.1) );
	     b = 0;
	  }
	  else if (col >=0.5 && col < 0.6)
	  {
	     r = 0;
	     g = int(175.0 + 80.0 * ((col - 0.5) / 0.1) );
	     b = int(160.0 * ((col- 0.5) / 0.1));
	  }
	  else if (col >=0.6 && col < 0.7)
	  {
	     r = int(120.0 * ((col - 0.6) / 0.1));
	     g = 255;
	     b = int(160.0 + 95.0 * ((col - 0.6) / 0.1) );
	  }
	  else if (col >=0.7 && col < 0.8)
	  {
	     r = int(120.0 - 50.0 * ((col - 0.7) / 0.1) );
	     g = int(255.0 - 185.0 * ((col - 0.7) / 0.1) );
	     b = 255;
	  }
	  else if (col >=0.8 && col < 0.9)
	  {
	     r = int(70.0 + 120.0 * ((col - 0.8) / 0.1) );
	     g = 70;
	     b = 255;
	  }
	  else if (col >=0.9 && col <= 1)
	  {
	     r = int(190 + 65.0 * ((col - 0.9) / 0.1) );
	     g = int(70.0 - 70.0 * ((col - 0.9) / 0.1) );
	     b = 255;
	  }
	}



	//an arctan function with outputs -pi/2 to 3pi/2
	float etan(float e_x, float e_y)
	{
	   float e_theta;
	   
	   if (e_x > 0)
	      e_theta = atan(e_y / e_x);
	   else if (e_x < 0)
	      e_theta = atan(e_y / e_x) + pi;
	   else if (e_x == 0 && e_y  > 0)
	      e_theta = pi / 2;
	   else if (e_x == 0 && e_y  < 0)
	      e_theta = -pi / 2;
	   else if (e_x == 0 && e_y == 0)
	      e_theta = 0;
	   
	   return e_theta;
	}



	//true if its closer to the plane
	bool is_it_closer(int x_co, int y_co, float distance)
	{
	   //is it closer constant default false
	   bool i_c = false;
	   
	   //if there is no point there
	   if (distance_to_the_plane[x_co][y_co] == -1)
	   {
	      distance_to_the_plane[x_co][y_co] = distance;
	      i_c = true;
	   }
	   //if the point is closer
	   if (distance_to_the_plane[x_co][y_co] != 0)
	   {
	      if (distance > distance_to_the_plane[x_co][y_co])
	      {
	         distance_to_the_plane[x_co][y_co] = distance;
	         i_c = true;
	      }
	   }
	   
	   return i_c;
	}



	//outputs an x coordinated based on the x, y, z coordinates and the rotations
	int determine_x_coord(float x_var, float y_var, float z_var,
	                      float theta_rotation, float phi_rotation, float sc_rotation)
	{
	   float t_rad, t_theta;
	   float x_c, t_x, t_x2;
	   float t_x1, t_y1;
	   
	   //theta rotation
	   t_rad = sqrt(pow(x_var,2) + pow(y_var,2));
	   t_theta = etan(x_var, y_var);
	   
	   t_theta-=theta_rotation;
	   
	   t_x1 = t_rad*cos(t_theta);
	   t_y1 = t_rad*sin(t_theta);
	   
	   //phi rotation
	   t_rad = sqrt(pow(t_x1,2) + pow(z_var,2));
	   t_theta = etan(t_x1, z_var);
	   
	   t_theta+=phi_rotation;
	   
	   t_x = t_rad * cos(t_theta);
	   
	   //screen rotation
	   t_rad =  sqrt(pow(t_x,2) + pow(t_y1,2));
	   t_theta = etan(t_x, t_y1);
	   
	   t_theta += sc_rotation; 
	   
	   t_x2 = t_rad*cos(t_theta);

	   x_c = (SCREEN_W / 2) + t_x2;
	    
	   return int(x_c);
	}



	//outputs a y coordinated based on the x, y, z coordinates and the rotations
	int determine_y_coord(float x_var, float y_var, float z_var,
	                      float theta_rotation, float phi_rotation, float sc_rotation)
	{  
	   float t_rad, t_theta;
	   float y_c, t_y, t_x;
	   float t_x1, t_y1;
	   
	   //theta rotation
	   t_rad = sqrt(pow(x_var,2) + pow(y_var,2));
	   t_theta = etan(x_var, y_var);
	   
	   t_theta-=theta_rotation;
	   
	   t_x1 = t_rad*cos(t_theta);
	   t_y1 = t_rad*sin(t_theta);
	   
	   //phi rotation
	   t_rad = sqrt(pow(t_x1,2) + pow(z_var,2));
	   t_theta = etan(t_x1, z_var);
	   
	   t_theta+=phi_rotation;
	   
	   t_x = t_rad * cos(t_theta);
	   
	   //screen rotation
	   t_rad =  sqrt(pow(t_x,2) + pow(t_y1,2));
	   t_theta = etan(t_x, t_y1);
	   
	   t_theta += sc_rotation;
	   
	   t_y = t_rad*sin(t_theta); 

	   y_c =(SCREEN_H / 2) -  t_y;
	   
	   return int(y_c);
	}



	//draws the axes based on the rotations
	void draw_axes(float theta_rotation, float phi_rotation, float scr_rotation)
	{
	   BITMAP *axes;
	   axes = create_bitmap(SCREEN_W, SCREEN_H);
	   clear_bitmap(axes);
	   
	   float p1, p2, p3, p4, p_norm, p_dist;
	   int x_coo, y_coo;
	   float a_x, a_y, a_z;
	   
	   p1 = sin(phi_rotation)*cos(theta_rotation);
	   p2 = sin(phi_rotation)*sin(theta_rotation);                    
	   p3 = cos(phi_rotation);
	   p4 = -1920;
	   p_norm = sqrt((p1 * p1) + (p2 * p2) + (p3 * p3));

	   
	   //x
	   for (a_x = -600; a_x < 600; a_x += 1)
	   {

	      x_coo = determine_x_coord(a_x, 0, 0, theta_rotation, phi_rotation, scr_rotation);
	      y_coo = determine_y_coord(a_x, 0, 0, theta_rotation, phi_rotation, scr_rotation);
	            
	      p_dist = ((p1 * a_x) + (p2 * 0) + (p3 * 0) + p4) / p_norm;
	      
	      if ((x_coo >= 0) && (x_coo <=1920) && 
	          (y_coo >=0) && (y_coo <=1200) && 
	          (is_it_closer(x_coo, y_coo, p_dist)) == true)
	      {
	         putpixel(axes, x_coo, y_coo, makecol(255, 0, 0));
	      } 
	   }
	   
	   //y
	   for (a_y = -600; a_y < 600; a_y += 1)
	   {

	      x_coo = determine_x_coord(0, a_y, 0, theta_rotation, phi_rotation, scr_rotation);
	      y_coo = determine_y_coord(0, a_y, 0, theta_rotation, phi_rotation, scr_rotation);
	  
	      p_dist = ((p1 * 0) + (p2 * a_y) + (p3 * 0) + p4) / p_norm;          
	   
	      if ((x_coo >= 0) && (x_coo <=1920) && 
	          (y_coo >=0) && (y_coo <=1200) && 
	          (is_it_closer(x_coo, y_coo, p_dist)) == true)
	      {
	         putpixel(axes, x_coo, y_coo, makecol(0, 255, 0));
	      } 
	   }
	   
	   //z
	   for (a_z = -600; a_z < 600; a_z += 1)
	   {

	      x_coo = determine_x_coord(0, 0, a_z, theta_rotation, phi_rotation, scr_rotation);
	      y_coo = determine_y_coord(0, 0, a_z, theta_rotation, phi_rotation, scr_rotation);
	            
	      p_dist = ((p1 * 0) + (p2 * 0) + (p3 * a_z) + p4) / p_norm;
	   
	      if ((x_coo >= 0) && (x_coo <=1920) && 
	          (y_coo >=0) && (y_coo <=1200) && 
	          (is_it_closer(x_coo, y_coo, p_dist)) == true)
	      {
	         putpixel(axes, x_coo, y_coo, makecol(0, 0, 255));
	      } 
	   }
	   
	   masked_blit(axes,buffer, 0, 0, 0, 0, SCREEN_W, SCREEN_H);
	   
	   destroy_bitmap(axes);
	}



	//draws the graph based on rotations coordinate system detail and zoom
	void draw_graph(int what_system, float t_rotation, float p_rotation, float s_rotation, float s_step, float scale, int which_one)
	{
	   //color
	   int r_c, g_c, b_c;  
	     
	   //percent complete  
	   int l_perc=0;  
	   
	   //a plane in the form (c1)x + (c2)y + (c3)z + (c4) = 0  
	   float c1, c2, c3, c4;
	   
	   //the norm of the vector <c1, c2, c3>, c_dist the distance to the plane   
	   float c_norm, c_dist; 
	   
	   
	      //spherical variables
	   float v_phi, v_theta, v_rho;
	   
	   //rectangular variables and polar
	   float v_x, v_y, v_z, v_r;
	   
	   
	   //2d coordinates
	   int x_coord, y_coord;
	   
	   //creates a plane based on theta and phi rotations
	   //sin(phi)cos(theta)x + sin(phi)sin(theta)y + cos(phi)z - rho = 0
	   c1 = sin(p_rotation)*cos(t_rotation);
	   c2 = sin(p_rotation)*sin(t_rotation);                    
	   c3 = cos(p_rotation);
	   c4 = -1920;
	   
	   //a vector normal to the plane
	   c_norm = sqrt((c1 * c1) + (c2 * c2) + (c3 * c3));
	   
	   //set the default distance to the plane to -1
	   for (int d1 = 0; d1 < 1920; d1++)
	   {
	      for (int d2 = 0; d2 < 1200; d2++)
	         distance_to_the_plane[d1][d2] = -1;                   
	       
	   }//endfor

	   if (s_step<0.02) s_step = s_step / 2;             

	   //phi from 0 to pi
	   for (v_phi = 0; v_phi < pi ; v_phi+=(s_step/2))
	   {
	      //theta from 0 to two pi
	      for (v_theta = 0; v_theta < 2*pi; v_theta+=s_step)
	      {
	          
	         //spherical coordinate system 
	         if (what_system == 1)
	         {
	            
	            //find rho as a function of theta and phi
	            if(which_one == 1)
	            {
	               //find rho
	               v_rho=what_rho_is(v_phi, v_theta, which_one, scale);
	            
	               //spherical to rectangular
	               v_x = scale * v_rho*sin(v_phi)*cos(v_theta);
	               v_y = scale * v_rho*sin(v_phi)*sin(v_theta);
	               v_z = scale * v_rho*cos(v_phi);
	               
	            }//endif
	            
	            //find theta as a function of phi and rho
	            if(which_one == 2)
	            {
	               //find theta (rho = theta, theta = rho)
	               v_rho=what_rho_is(v_phi, v_theta, which_one, scale);
	            
	               //spherical to rectangular
	               v_x = (v_theta / (2 * pi)) * (750)*sin(v_phi)*cos(v_rho);
	               v_y = (v_theta / (2 * pi)) * (750)*sin(v_phi)*sin(v_rho);
	               v_z = (v_theta / (2 * pi)) * (750)*cos(v_phi);
	               
	            }//endif
	            
	            //find phi as a function of theta and rho
	            if(which_one == 3)
	            {
	               //find phi (rho = phi, theta = phi)
	               v_rho=what_rho_is(v_phi, v_theta, which_one, scale);
	            
	               //spherical to rectangular
	               v_x = (v_phi / pi) * (750) * sin(v_rho)*cos(v_theta);
	               v_y = (v_phi / pi) * (750) * sin(v_rho)*sin(v_theta);
	               v_z = (v_phi / pi) * (750) * cos(v_rho);
	               
	            }//endif
	            
	            
	               
	         }//endif spherical
	          
	            
	         //rectangular coordinates
	         if (what_system == 2)
	         { 
	            //z = f(x, y)
	            if(which_one == 1)
	            {
	               //change parameters
	               //0 -> 2pi to -1000 -> 1000
	               v_x = (v_theta * 750 / pi - 750);
	               //0 -> pi to -1000 -> 1000
	               v_y = (v_phi * 1500 / pi - 750);
	               
	               //find z
	               v_z = scale * what_z_is(v_x/scale, v_y/scale, which_one);
	               
	            }//endif
	            
	            //y = f(x,z)
	            if(which_one == 2)
	            {
	               //change parameters
	               //0 -> 2pi to -1000 -> 1000
	               v_x = (v_theta * 1000 / pi - 1000);
	               //0 -> pi to -1000 -> 1000
	               v_z = (v_phi * 2000 / pi - 1000);
	               
	               //find z
	               v_y = scale * what_z_is(v_x/scale, v_z/scale, which_one);
	               
	            }//endif
	            
	            //x = f(z,y)
	            if(which_one == 3)
	            {
	               //change parameters
	               //0 -> 2pi to -1000 -> 1000
	               v_z = (v_theta * 1000 / pi - 1000);
	               //0 -> pi to -1000 -> 1000
	               v_y = (v_phi * 2000 / pi - 1000);
	               
	               //find z
	               v_x = scale * what_z_is(v_z/scale, v_y/scale, which_one);
	               
	            }//endif
	               
	            
	         }//endif rectangular 
	         
	         
	         //parametric
	         if (what_system == 3)
	         { 
	            //if(which_one == 1)
	            //{
	               //all 3 at once with pointers
	               parametric_xyz(v_x, v_y, v_z, v_phi, v_theta);
	               
	               //blow up
	               v_x = scale * v_x;
	               v_y = scale * v_y;
	               v_z = scale * v_z;
	               
	            //}//endif
	               
	            
	         }//endif
	        
	         
	         //cylindrical
	         if (what_system == 4)
	         { 
	            if(which_one == 1)
	            {
	             
	               v_r = (v_phi * 500) / pi ;      
	                       
	               //find z
	               v_z = scale * what_polar_z_is(v_r, v_theta, which_one,scale)/scale;

	               
	               //find x, y
	               v_x = scale *(v_r * cos(v_theta));
	               v_y = scale *(v_r * sin(v_theta));
	               
	            }//endif
	            
	            if(which_one == 2)
	            {
	              //phi = r 0 -> pi to 0 -> 750 
	              v_z = (v_phi/ pi) * (750) -(750/2);
	               
	              //find r 
	               v_r = scale * what_polar_z_is(v_z, v_theta, which_one,scale)/scale;
	               
	               
	               
	               //find x, y
	               v_x = scale * (v_r * cos(v_theta));
	               v_y = scale * (v_r * sin(v_theta));
	               
	            }//endif
	            
	            if(which_one == 3)
	            {
	              v_z = (v_theta/(2* pi)) * (750) -(750/2);
	              v_rho = (v_phi * 500) / pi ;     
	                
	              //find theta v_r = theta 
	               v_r = scale * what_polar_z_is(v_rho, v_z , which_one,scale)/scale;
	               
	              
	               //find x, y
	               v_x = scale * (v_rho) * cos(v_r);
	               v_y = scale * (v_rho) * sin(v_r);
	               
	            }//endif
	               
	            
	         }//endif
	            
	         if (v_z < 1920 && v_x < 1920 && v_y < 1920)
	         {
	            //find the 2d projection
	            x_coord = determine_x_coord(v_x, v_y, v_z, t_rotation, p_rotation, s_rotation);
	            y_coord = determine_y_coord(v_x, v_y, v_z, t_rotation, p_rotation, s_rotation);
	         
	            //find the distance to the plane
	            c_dist = ((c1 * v_x) + (c2 * v_y) + (c3 * v_z) + c4) / c_norm;
	         
	           e_color(r_c, g_c, b_c, 1-sin(v_phi) );
	      
	            if ((x_coord >= 0) && (x_coord <=1920) && 
	                (y_coord >=0) && (y_coord <=1200) && 
	                (is_it_closer(x_coord, y_coord, c_dist)) == true)
	            {
	               putpixel(buffer, x_coord, y_coord, makecol(r_c, g_c, b_c));
	            } 
	              
	         } //end if (v_z)
	            
	      } //end for (theta)
	      l_perc=int((100*v_phi)/pi);   
	      textprintf_ex(screen, font, 20, 20, makecol(255, 255, 255),  bitmap_mask_color(screen), "loading %-5d...", l_perc);         
	   } //end for (phi)
	   
	}



	//used to draw buttons for coordinate systems
	void draw_coord_buttons(int what_sys, int o_sys)
	{
	     
	   //redraw circles for what system sph rect para
	   circle(buffer, 1500, 720, 15, makecol(255,255,255));
	   circle(buffer, 1680, 720, 15, makecol(255,255,255));
	   circle(buffer, 1850, 720, 15, makecol(255,255,255));
	   circle(buffer, 1590, 720, 15, makecol(255,255,255));
	   //circle(buffer, 1770, 720, 15, makecol(255,255,255));
	     
	   if (o_sys == 1)
	      circlefill(buffer, 1500, 720, 10, bitmap_mask_color(screen));
	   if (o_sys == 2)   
	      circlefill(buffer, 1680, 720, 10, bitmap_mask_color(screen));
	   if (o_sys == 3)   
	      circlefill(buffer, 1850, 720, 10, bitmap_mask_color(screen));
	   if (o_sys == 4)   
	      circlefill(buffer, 1590, 720, 10, bitmap_mask_color(screen));
	  // if (o_sys == 5)   
	  //    circlefill(buffer, 1770, 720, 10, bitmap_mask_color(screen));
	   
	   if (what_sys == 1)
	      circlefill(buffer, 1500, 720, 10, makecol(255,255,255));
	   if (what_sys == 2)
	      circlefill(buffer, 1680, 720, 10,  makecol(255,255,255));
	   if (what_sys == 3)
	      circlefill(buffer, 1850, 720, 10, makecol(255,255,255));
	   if (what_sys == 4)   
	      circlefill(buffer, 1590, 720, 10,  makecol(255,255,255));
	  // if (what_sys == 5)   
	  //    circlefill(buffer, 1770, 720, 10,  makecol(255,255,255));
	     
	}



	///draws the interface based on the positions of the buttons
	void draw_interface(int mx, int mx2, int mx3, int mx4, int mx5, int w_sys, int o_w_sys, bool dy)
	{    
	   //redraw bars detail/zoom/phi/theta/screen/redraw/rotate
	   masked_blit (bar, buffer,0,0, 1480,780,1880,820);
	   masked_blit (bar, buffer,0,0, 1480,840,1880,880);
	   masked_blit (bar, buffer,0,0, 1480,900,1880,940);
	   masked_blit (bar, buffer,0,0, 1480,960,1880,1000); 
	   masked_blit (bar, buffer,0,0, 1480,1020,1880,1060);
	   masked_blit (rota2, buffer,0,0, 1730,1100,1830,1160);    
	                    
	   //change button color
	   if (dy)
	   {
	      blit (rdra2, buffer,0,0, 1530,1100,1630,1160);
	      dy = false;
	   }
	   else  
	   {       
	      blit (rdra, buffer,0,0, 1530,1100,1630,1160);
	      dy = true;
	   }

	   draw_coord_buttons(w_sys, o_w_sys);
	               
	   //redraw boxes on bars
	   blit(box,buffer, 0,0, mx4-8,784,mx4+8,816);
	   blit(box,buffer, 0,0, mx5-8,844,mx5+8,876); 
	   blit(box,buffer, 0,0, mx3-8,904,mx3+8,936);
	   blit(box,buffer, 0,0, mx-8,964,mx+8,996); 
	   blit(box,buffer, 0,0, mx2-8,1040-16,mx2+8,1040+16);    
	}



	//w_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash)
	void draw_everything(int wh_s, int owh_s, float t_ro, float p_ro, float sc_ro, float de, float zo, int tr_b, int pr_b, int sr_b, int d_b, int s_b, bool b_fl, int w_t)
	{
	   //draw everything
	   clear_bitmap(buffer);
	   draw_graph(wh_s, t_ro, p_ro, sc_ro, de, zo, w_t);
	   draw_interface(tr_b, sr_b, pr_b, d_b, s_b, wh_s, owh_s, b_fl);
	   draw_axes(t_ro, p_ro, sc_ro);
	   textprintf_ex(buffer, font, 1490, 675, makecol(255, 255, 255),  bitmap_mask_color(screen), "sph        pol        rect                  par  ");
	   textprintf_ex(buffer, font, 1490, 765, makecol(255, 255, 255),  bitmap_mask_color(screen), "less          detail =  %-5f            more", de);
	   textprintf_ex(buffer, font, 1490, 825, makecol(255, 255, 255),  bitmap_mask_color(screen), "less        zoom    =    %-5f          more", 5 * pow(10, ((float(s_b - 1500) / 90) - 2) ));
	   textprintf_ex(buffer, font, 1490, 885, makecol(255, 255, 255),  bitmap_mask_color(screen), "0             phi rotation =  %-5d           pi", (pr_b-1500)/2);
	   textprintf_ex(buffer, font, 1490, 945, makecol(255, 255, 255),  bitmap_mask_color(screen), "0            theta rotation = %-5d          2pi", tr_b-1500);
	   textprintf_ex(buffer, font, 1490, 1005, makecol(255, 255, 255),  bitmap_mask_color(screen),"0           screen rotation = %-5d          2pi", sr_b-1500);               
	   blit(buffer, screen, 0, 0, 0, 0, SCREEN_W, SCREEN_H);
	   vsync(); 
	}
	   

	int main()
	{ 
	  //initialize allegro
	   allegro_init();
	   install_keyboard();  
	   install_mouse();    
	                                                                   
	   if (set_gfx_mode(GFX_AUTODETECT, 1920, 1200, 0, 0))                                      
	   {                                                      
	      allegro_message("Error setting GFX mode");                                             
	      return 1;   
	   }                                     
	   
	   
	   //rotations size and scale
	   float detail = pi/90; 
	   float zoom = 400;
	   float th_rotation=pi/4, ph_rotation=pi/4, sc_rotation=3*pi/2;
	   //for button color switch
	   bool button_flash = true;   
	   //for displaying rotations ect
	   int t_rot_bar_x=1545,s_rot_bar_x=1770,p_rot_bar_x=1590,d_bar_x=1500,s_bar_x=1700;
	   //what system 1 sph 2 rect 3 para and proir used system
	   int w_s=1, ol_s=1, which_t=1;
	   
	   
	   
	   //palette
	   //set_palette(desktop_palette);
	   
	   //create a bitmap for double buffering equal to the screen size declared globally
	   buffer = create_bitmap(SCREEN_W, SCREEN_H);
	   clear_bitmap(buffer);
	   //buffer2 used for the mouse
	   BITMAP *buffer2;
	   buffer2 = create_bitmap(SCREEN_W, SCREEN_H);
	   clear_bitmap(buffer2);
	   
	   
	   
	   //bitmap used for mouse cursor
	   BITMAP *mouse_cursor;
	   mouse_cursor = create_bitmap(27,27);
	   clear_to_color(mouse_cursor, bitmap_mask_color(screen));
	   for (float c=0; c<10; c+=.005)
	      circle(mouse_cursor, 13, 13, int(c/2), makecol(255-int(255-16*c),0,0));
	   set_mouse_sprite(mouse_cursor);
	   set_mouse_sprite_focus(13, 13);
	   


	   //bar
	   bar = create_bitmap(400, 40);
	   clear_bitmap(bar);
	   for (int k2 =20;k2 < 200; k2++)
	   {
	      circle(bar, 400-k2, 20, 5, makecol(int((sin((k2-20)*3.15159/180)+1)*63),int((sin((k2-20)*3.15159/180)+1)*255),int((sin((k2-20)*3.15159/1800)+1)*255)));
	      circle(bar, k2, 20, 5, makecol(int((sin((k2-20)*3.15159/180)+1)*63),int((sin((k2-20)*3.15159/180)+1)*255),int((sin((k2-20)*3.15159/1800)+1)*255)));
	   }
	   line(bar,20,20,380,20,makecol(0,200,200));
	   
	         
	   //nob
	   box = create_bitmap(16,32);
	   rectfill(box, 0,0,16,32,makecol(255,63,127));
	   rectfill(box, 3,3,12,28,makecol(55,255,255));
	   
	   
	   
	   //redraw 1
	   rdra = create_bitmap(100, 60);
	   //redraw 2
	   rdra2 = create_bitmap(100, 60);
	   //rotate 1
	   rota = create_bitmap(100, 60);
	   //rotate 2
	   rota2 = create_bitmap(100, 60);
	   
	   
	   
	   //draw 2x2 buttons
	   for(int k3=0; k3<20; k3++)
	   {
	      rect( rdra, k3, k3, 100 - k3, 60 - k3, makecol(55 - k3, k3 * 5, k3 + 75) );
	      textprintf_ex(rdra, font, 27, 27, makecol(123, 123, 231), bitmap_mask_color(screen), "REDRAW");
	      rect( rdra2, k3, k3, 100 - k3, 60 - k3, makecol(85 - k3,k3,5 * k3) );
	      textprintf_ex(rdra2, font, 27, 27, makecol(123, 231, 123), bitmap_mask_color(screen), "REDRAW");
	      rect( rota, k3, k3, 100 - k3, 60 - k3, makecol(55 - k3, k3 * 5, k3 + 75) );
	      textprintf_ex(rota, font, 27, 27, makecol(123, 231,123), bitmap_mask_color(screen), " SAVE ");
	      rect( rota2, k3, k3, 100 - k3, 60 - k3, makecol(85 - k3, k3, 5 * k3) );
	      textprintf_ex(rota2, font, 27, 27, makecol(123, 231,123), bitmap_mask_color(screen), " SAVE ");
	   }
	     
	     
	   draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
	   
	   do
	   {
	      //left click
	      if (mouse_b & 1)
	      {
	         
	         
	         //if the detail bar is clicked on
	         
	         if ((mouse_x>= 1500) && (mouse_x <= 1860) && (mouse_y >= 760) && (mouse_y <= 820))
	         {
	            detail=mouse_x;
	            d_bar_x=mouse_x;
	            
	            blit (bar, buffer,0,0, 1480,780,1880,820);
	            blit(box, buffer, 0,0, d_bar_x-8,784,d_bar_x+8,816);
	            
	            //0 to pi/720
	            detail=( (-7*pi/720) * (  (detail-1500) /360  ) )+(pi/90);
	            
	            //for fast rotation
	            if (detail > 0.02)
	               draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);  
	         }   
	         
	         
	         //if the size bar is clicked on
	         
	         if ((mouse_x>= 1500) && (mouse_x <= 1860) && (mouse_y >= 840) && (mouse_y <= 880))
	         {
	            zoom=mouse_x;
	            s_bar_x=mouse_x;
	            
	            blit (bar, buffer,0,0, 1480,840,1880,880);
	            blit(box,buffer, 0,0, s_bar_x-8,844,s_bar_x+8,876); 
	            
	            //0 to 720
	            zoom = 5 * pow(10, ((float(s_bar_x - 1500) / 90) - 2) );
	            
	            //for fast rotation
	            if (detail > 0.02)
	               draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t); 
	         } 
	                 
	                 
	         //if the phi rot bar is clicked on        
	         
	         if ((mouse_x>= 1500) && (mouse_x <= 1860) && (mouse_y >= 900) && (mouse_y <= 940))
	         {
	            ph_rotation=mouse_x;
	            p_rot_bar_x=mouse_x;

	            blit (bar, buffer,0,0, 1480,900,1880,940);
	            blit(box,buffer, 0,0, p_rot_bar_x-8,904,p_rot_bar_x+8,936); 
	            
	            //0 to pi
	            ph_rotation = ((ph_rotation-1500)*pi)/360;
	            
	            //for fast rotation
	            if (detail > 0.02)
	               draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t); 
	         }
	         
	         
	         //if the theta rot bar is clicked on
	         
	         if ((mouse_x>= 1500) && (mouse_x <= 1860) && (mouse_y >= 960) && (mouse_y <= 1000))
	         {
	            th_rotation=mouse_x;
	            t_rot_bar_x=mouse_x;

	            blit (bar, buffer,0,0, 1480,960,1880,1000);
	            blit(box,buffer, 0,0, t_rot_bar_x-8,964,t_rot_bar_x+8,996); 
	            
	            //0 to 2pi
	            th_rotation = ((th_rotation-1500)*pi)/180;
	            
	            //for fast rotation
	            if (detail > 0.02)
	               draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t); 
	         }
	         
	         // if screen rotation is clicked
	         
	         
	         if ((mouse_x>= 1500) && (mouse_x <= 1860) && (mouse_y >= 1020) && (mouse_y <= 1060))
	         {
	            sc_rotation=mouse_x;
	            s_rot_bar_x=mouse_x;            
	            
	            blit (bar, buffer,0,0, 1480,1020,1880,1060);
	            blit(box,buffer, 0,0, s_rot_bar_x-8,1040-16,s_rot_bar_x+8,1040+16); 
	            
	            //0 to 2pi
	            sc_rotation = ((sc_rotation-1500)*pi)/180;
	            
	            //for fast rotation
	            if (detail > 0.02)
	               draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t); 
	         }
	         
	         
	         
	         //redraw
	         if ((mouse_x>= 1530) && (mouse_x <= 1630) && (mouse_y >= 1110) && (mouse_y <= 1160))
	         {
	            //flash button
	            if (button_flash)
	               button_flash = false;
	            else        
	               button_flash = true;
	            //draw graph interface etc
	            draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);          
	         }
	         
	         
	         
	         //go / rotate
	         if ((mouse_x>= 1730) && (mouse_x <= 1830) && (mouse_y >= 1110) && (mouse_y <= 1160))
	         {
	            save_bitmap("graph.bmp", buffer, 0);
	            if (button_flash)
	            {
	               blit (rota, buffer, 0, 0, 1730,1100,1830,1160);
	               button_flash = false;
	            }
	            else  
	            {       
	               blit (rota2, buffer,0,0, 1730,1100,1830,1160);
	               button_flash = true;
	            }
	         }
	         
	         
	         //spherical
	         if (pow(mouse_x-1500,2)+pow(mouse_y-720,2)<=225)
	         {
	            ol_s=w_s;
	            w_s=1;
	            if (detail > 0.02)
	               draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t); 
	            else   
	               draw_coord_buttons(w_s, ol_s);
	         }
	         //rectangular
	         if (pow(mouse_x-1680,2)+pow(mouse_y-720,2)<=225)
	         {
	            ol_s=w_s;
	            w_s = 2;
	            if (detail > 0.02)
	               draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t); 
	            else   
	               draw_coord_buttons(w_s, ol_s);
	         }
	         //parametric
	         if (pow(mouse_x-1850,2)+pow(mouse_y-720,2)<=225)
	         {
	            ol_s=w_s;
	            w_s=3;
	            if (detail > 0.02)
	               draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t); 
	            else   
	               draw_coord_buttons(w_s, ol_s);
	         }  
	          //polar
	         if (pow(mouse_x-1590,2)+pow(mouse_y-720,2)<=225)
	         {
	            ol_s=w_s;
	            w_s=4;
	            if (detail > 0.02)
	               draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t); 
	            else   
	               draw_coord_buttons(w_s, ol_s);
	         }  
	          //diff eq
	         if (pow(mouse_x-1770,2)+pow(mouse_y-720,2)<=225)
	         {
	            ol_s=w_s;
	            w_s=5;
	            if (detail > 0.02)
	               draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t); 
	            else   
	               draw_coord_buttons(w_s, ol_s);
	         }  
	         
	      }
	      
	      //right click
	      if (mouse_b & 2)
	      {
			 which_t++;
			 if (which_t >=4)
			    which_t = 1;
			 draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
	      }

	      clear_bitmap(buffer2);
	      show_mouse(buffer2);
	      

	      blit(buffer, screen, 0, 0, 0, 0, SCREEN_W, SCREEN_H);
	      masked_blit(buffer2, screen, 0, 0, 0, 0, SCREEN_W, SCREEN_H);
	      vsync();
	            
	   }while(!key[KEY_ESC]);

	  
	  
	  
	   destroy_bitmap(buffer);
	   destroy_bitmap(buffer2);
	   destroy_bitmap(mouse_cursor);
	   destroy_bitmap(bar);
	   destroy_bitmap(box);
	   destroy_bitmap(rdra);
	   destroy_bitmap(rdra2);
	   
	   return 0;
	   
	}END_OF_MAIN();


}
