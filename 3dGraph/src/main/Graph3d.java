package main;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JPanel;

public class Graph3d {
	private static final int SCREEN_W = 1920;
	private static final int SCREEN_H = 1200;

	private static double pi = 3.14159265;
	private static double two_pi = 2 * pi;
	//array used by main and is_it_closer function
	private static double[][] distance_to_the_plane = new double[1921][1201];
	//what_system   1 = spherical 2 = rectangular 4 = polar
	//where the graph and interface are drawn
	private static BufferedImage bufferImage;
	private static Graphics buffer;
	private static BufferedImage mouse_cursorImage;

	//for interface
	private static BufferedImage bar;
	private static BufferedImage box;
	private static BufferedImage rdra;
	private static BufferedImage rdra2;
	private static BufferedImage rota;
	private static BufferedImage rota2;

	private static Graphics screen;
	private static Font font = new Font("consolas", Font.PLAIN, 14);

	//the equation for spherical coordinates
	double what_rho_is(double ph, double th, int which_type, double how_large) {
		double rho, theta, phi;
		double returned_var = 0.0;

		if (which_type == 1) {
			theta = th;
			phi = ph;
			/*-(rho)-*/
			rho = 1 - Math.sqrt(Math.tan(Math.sin(phi)));
			/*-(rho)-*/
			returned_var = rho;
		} else if (which_type == 2) {
			phi = ph;
			rho = (th / two_pi) * (750 / (how_large + 1));
			/*-(theta)-*/
			theta = phi - pi * Math.sin(rho);
			/*-(theta)-*/
			returned_var = theta;
		} else if (which_type == 3) {
			theta = th;
			rho = (ph / pi) * (750 / (how_large + 1));
			/*-(phi)-*/
			phi = 1 - rho * Math.sin(theta);
			/*-(phi)-*/
			returned_var = phi;
		}

		return returned_var;
	}

	//the equation for rectangular coordinates
	double what_z_is(double x1, double y1, int which_type) {
		double x, y, z;
		double returned_var = 0.0;

		if (which_type == 1) {
			x = x1;
			y = y1;
			/*-(z)-*/
			z = Math.sin(Math.sqrt(x * y));
			/*-(z)-*/
			returned_var = z;
		} else if (which_type == 2) {
			x = x1;
			z = y1;
			/*-(y)-*/
			y = x * x * x - z * z;
			/*-(y)-*/
			returned_var = y;
		} else if (which_type == 3) {
			z = x1;
			y = y1;
			/*-(x)-*/
			x = z + y * y - Math.sin(z);
			/*-(x)-*/
			returned_var = x;
		}

		return returned_var;
	}

	//the equation in polar
	double what_polar_z_is(double phi_va, double theta_va, int which_type, double how_large) {
		double z, r, theta;
		double returned_var = 0.0;

		if (which_type == 1) {
			theta = theta_va;
			r = phi_va;
			/*-(z)-*/
			z = Math.sqrt(r);
			/*-(z)-*/
			returned_var = z;
		} else if (which_type == 2) {
			theta = theta_va;
			z = phi_va; //(phi_va / pi) * (750/(how_large)) -((750/(how_large))/2);
			/*-(r)-*/
			r = Math.sin(z - theta);
			/*-(r)-*/
			returned_var = r;
		} else if (which_type == 3) {
			r = phi_va;//(phi_va/ pi ) * (1500/how_large);
			z = theta_va; //(theta_va / two_pi) * (750/(how_large)) -((750/(how_large))/2);
			/*-(r)-*/
			theta = z - Math.sin(Math.cos(1 - r));
			/*-(r)-*/
			returned_var = theta;
		}

		return returned_var;
	}

	private static class ParametricResult {
		final double x_uv;
		final double y_uv;
		final double z_uv;

		public ParametricResult(double x_uv, double y_uv, double z_uv) {
			this.x_uv = x_uv;
			this.y_uv = y_uv;
			this.z_uv = z_uv;
		}
	}

	//parametric equations
	ParametricResult parametric_xyz(double phi_u, double theta_v) {
		double u_min = 0;
		double u_max = 2 * pi;
		double v_min = 0;
		double v_max = 2 * pi;
		double u, v;

		//change parameters
		u = (phi_u / pi) * (u_max - u_min) + u_min;
		v = (theta_v / two_pi) * (v_max - v_min) + v_min;

		/*-(parametric)-*/
		double x_uv = Math.cos(v);
		double y_uv = Math.sin(v) * Math.cos(u);
		double z_uv = Math.sin(u);
		/*-(parametric)-*/

		return new ParametricResult(x_uv, y_uv, z_uv);
	}

	Color e_color(double col) {
		double r = 0.0;
		double g = 0.0;
		double b = 0.0;

		if (col >= 0 && col < 0.1) {
			r = 255;
			g = 145.0 * (col / 0.1);
			b = 30.0 * (col / 0.1);
		} else if (col >= 0.1 && col < 0.2) {
			r = 255;
			g = 145.0 + 110.0 * ((col - 0.1) / 0.1);
			b = 30.0 + 20.0 * ((col - 0.1) / 0.1);
		} else if (col >= 0.2 && col < 0.3) {
			r = 255.0 - 80.0 * ((col - 0.2) / 0.1);
			g = 255;
			b = 50.0 - 50.0 * ((col - 0.2) / 0.1);
		} else if (col >= 0.3 && col < 0.4) {
			r = 175.0 - 175.0 * ((col - 0.3) / 0.1);
			g = 255;
			b = 0;
		} else if (col >= 0.4 && col < 0.5) {
			r = 0;
			g = 255.0 - 80.0 * ((col - 0.4) / 0.1);
			b = 0;
		} else if (col >= 0.5 && col < 0.6) {
			r = 0;
			g = 175.0 + 80.0 * ((col - 0.5) / 0.1);
			b = 160.0 * ((col - 0.5) / 0.1);
		} else if (col >= 0.6 && col < 0.7) {
			r = 120.0 * ((col - 0.6) / 0.1);
			g = 255;
			b = 160.0 + 95.0 * ((col - 0.6) / 0.1);
		} else if (col >= 0.7 && col < 0.8) {
			r = 120.0 - 50.0 * ((col - 0.7) / 0.1);
			g = 255.0 - 185.0 * ((col - 0.7) / 0.1);
			b = 255;
		} else if (col >= 0.8 && col < 0.9) {
			r = 70.0 + 120.0 * ((col - 0.8) / 0.1);
			g = 70;
			b = 255;
		} else if (col >= 0.9 && col <= 1) {
			r = 190 + 65.0 * ((col - 0.9) / 0.1);
			g = 70.0 - 70.0 * ((col - 0.9) / 0.1);
			b = 255;
		}

		return makecol((int) r, (int) g, (int) b);
	}

	//an arcMath.tan function with outputs -pi/2 to 3pi/2
	double etan(double e_x, double e_y) {
		double e_theta = 0.0;

		if (e_x > 0) {
			e_theta = Math.atan(e_y / e_x);
		} else if (e_x < 0) {
			e_theta = Math.atan(e_y / e_x) + pi;
		} else if (e_x == 0 && e_y > 0) {
			e_theta = pi / 2;
		} else if (e_x == 0 && e_y < 0) {
			e_theta = -pi / 2;
		} else if (e_x == 0 && e_y == 0) {
			e_theta = 0;
		}

		return e_theta;
	}

	//true if its closer to the plane
	boolean is_it_closer(int x_co, int y_co, double distance) {
		//is it closer consMath.tant default false
		boolean i_c = false;

		//if there is no point there
		if (distance_to_the_plane[x_co][y_co] == -1) {
			distance_to_the_plane[x_co][y_co] = distance;
			i_c = true;
		}
		//if the point is closer
		if (distance_to_the_plane[x_co][y_co] != 0) {
			if (distance > distance_to_the_plane[x_co][y_co]) {
				distance_to_the_plane[x_co][y_co] = distance;
				i_c = true;
			}
		}

		return i_c;
	}

	//outputs an x coordinated based on the x, y, z coordinates and the rotations
	int determine_x_coord(double x_var, double y_var, double z_var, double theta_rotation, double phi_rotation, double sc_rotation) {
		double t_rad, t_theta;
		double x_c, t_x, t_x2;
		double t_x1, t_y1;

		//theta rotation
		t_rad = Math.sqrt(Math.pow(x_var, 2) + Math.pow(y_var, 2));
		t_theta = etan(x_var, y_var);

		t_theta -= theta_rotation;

		t_x1 = t_rad * Math.cos(t_theta);
		t_y1 = t_rad * Math.sin(t_theta);

		//phi rotation
		t_rad = Math.sqrt(Math.pow(t_x1, 2) + Math.pow(z_var, 2));
		t_theta = etan(t_x1, z_var);

		t_theta += phi_rotation;

		t_x = t_rad * Math.cos(t_theta);

		//screen rotation
		t_rad = Math.sqrt(Math.pow(t_x, 2) + Math.pow(t_y1, 2));
		t_theta = etan(t_x, t_y1);

		t_theta += sc_rotation;

		t_x2 = t_rad * Math.cos(t_theta);

		x_c = (SCREEN_W / 2) + t_x2;

		return (int) x_c;
	}

	//outputs a y coordinated based on the x, y, z coordinates and the rotations
	int determine_y_coord(double x_var, double y_var, double z_var, double theta_rotation, double phi_rotation, double sc_rotation) {
		double t_rad, t_theta;
		double y_c, t_y, t_x;
		double t_x1, t_y1;

		//theta rotation
		t_rad = Math.sqrt(Math.pow(x_var, 2) + Math.pow(y_var, 2));
		t_theta = etan(x_var, y_var);

		t_theta -= theta_rotation;

		t_x1 = t_rad * Math.cos(t_theta);
		t_y1 = t_rad * Math.sin(t_theta);

		//phi rotation
		t_rad = Math.sqrt(Math.pow(t_x1, 2) + Math.pow(z_var, 2));
		t_theta = etan(t_x1, z_var);

		t_theta += phi_rotation;

		t_x = t_rad * Math.cos(t_theta);

		//screen rotation
		t_rad = Math.sqrt(Math.pow(t_x, 2) + Math.pow(t_y1, 2));
		t_theta = etan(t_x, t_y1);

		t_theta += sc_rotation;

		t_y = t_rad * Math.sin(t_theta);

		y_c = (SCREEN_H / 2) - t_y;

		return (int) y_c;
	}

	//draws the axes based on the rotations
	void draw_axes(double theta_rotation, double phi_rotation, double scr_rotation) {
		double p1, p2, p3, p4, p_norm, p_dist;
		int x_coo, y_coo;
		double a_x, a_y, a_z;

		p1 = Math.sin(phi_rotation) * Math.cos(theta_rotation);
		p2 = Math.sin(phi_rotation) * Math.sin(theta_rotation);
		p3 = Math.cos(phi_rotation);
		p4 = -1920;
		p_norm = Math.sqrt((p1 * p1) + (p2 * p2) + (p3 * p3));

		//x
		for (a_x = -600; a_x < 600; a_x += 1) {
			x_coo = determine_x_coord(a_x, 0, 0, theta_rotation, phi_rotation, scr_rotation);
			y_coo = determine_y_coord(a_x, 0, 0, theta_rotation, phi_rotation, scr_rotation);

			p_dist = ((p1 * a_x) + (p2 * 0) + (p3 * 0) + p4) / p_norm;

			if ((x_coo >= 0) && (x_coo <= 1920) && (y_coo >= 0) && (y_coo <= 1200) && (is_it_closer(x_coo, y_coo, p_dist)) == true) {
				buffer.setColor(Color.RED);
				buffer.drawLine(x_coo, y_coo, x_coo, y_coo);
			}
		}

		//y
		for (a_y = -600; a_y < 600; a_y += 1) {
			x_coo = determine_x_coord(0, a_y, 0, theta_rotation, phi_rotation, scr_rotation);
			y_coo = determine_y_coord(0, a_y, 0, theta_rotation, phi_rotation, scr_rotation);

			p_dist = ((p1 * 0) + (p2 * a_y) + (p3 * 0) + p4) / p_norm;

			if ((x_coo >= 0) && (x_coo <= 1920) && (y_coo >= 0) && (y_coo <= 1200) && (is_it_closer(x_coo, y_coo, p_dist)) == true) {
				buffer.setColor(Color.GREEN);
				buffer.drawLine(x_coo, y_coo, x_coo, y_coo);
			}
		}

		//z
		for (a_z = -600; a_z < 600; a_z += 1) {
			x_coo = determine_x_coord(0, 0, a_z, theta_rotation, phi_rotation, scr_rotation);
			y_coo = determine_y_coord(0, 0, a_z, theta_rotation, phi_rotation, scr_rotation);

			p_dist = ((p1 * 0) + (p2 * 0) + (p3 * a_z) + p4) / p_norm;

			if ((x_coo >= 0) && (x_coo <= 1920) && (y_coo >= 0) && (y_coo <= 1200) && (is_it_closer(x_coo, y_coo, p_dist)) == true) {
				buffer.setColor(Color.BLUE);
				buffer.drawLine(x_coo, y_coo, x_coo, y_coo);
			}
		}
	}

	//draws the graph based on rotations coordinate system detail and zoom
	void draw_graph(int what_system, double t_rotation, double p_rotation, double s_rotation, double s_step, double scale, int which_one) {
		//percent complete
		int l_perc = 0;

		//a plane in the form (c1)x + (c2)y + (c3)z + (c4) = 0
		double c1, c2, c3, c4;

		//the norm of the vector <c1, c2, c3>, c_dist the distance to the plane
		double c_norm, c_dist;

		//spherical variables
		double v_phi, v_theta, v_rho;

		//rectangular variables and polar
		double v_r;
		double v_x = 0;
		double v_y = 0;
		double v_z = 0;

		//2d coordinates
		int x_coord, y_coord;

		//creates a plane based on theta and phi rotations
		//sin(phi)cos(theta)x + sin(phi)sin(theta)y + cos(phi)z - rho = 0
		c1 = Math.sin(p_rotation) * Math.cos(t_rotation);
		c2 = Math.sin(p_rotation) * Math.sin(t_rotation);
		c3 = Math.cos(p_rotation);
		c4 = -1920;

		//a vector normal to the plane
		c_norm = Math.sqrt((c1 * c1) + (c2 * c2) + (c3 * c3));

		//set the default disMath.tance to the plane to -1
		for (int d1 = 0; d1 < 1920; d1++) {
			for (int d2 = 0; d2 < 1200; d2++) {
				distance_to_the_plane[d1][d2] = -1;
			}
		} //endfor

		if (s_step < 0.02) {
			s_step = s_step / 2;
		}

		//phi from 0 to pi
		for (v_phi = 0; v_phi < pi; v_phi += (s_step / 2)) {
			//theta from 0 to two pi
			for (v_theta = 0; v_theta < 2 * pi; v_theta += s_step) {
				//spherical coordinate system
				if (what_system == 1) {
					//find rho as a function of theta and phi
					if (which_one == 1) {
						//find rho
						v_rho = what_rho_is(v_phi, v_theta, which_one, scale);
						//spherical to rectangular
						v_x = scale * v_rho * Math.sin(v_phi) * Math.cos(v_theta);
						v_y = scale * v_rho * Math.sin(v_phi) * Math.sin(v_theta);
						v_z = scale * v_rho * Math.cos(v_phi);
					} else if (which_one == 2) {//find theta as a function of phi and rho
						//find theta (rho = theta, theta = rho)
						v_rho = what_rho_is(v_phi, v_theta, which_one, scale);
						//spherical to rectangular
						v_x = (v_theta / (2 * pi)) * (750) * Math.sin(v_phi) * Math.cos(v_rho);
						v_y = (v_theta / (2 * pi)) * (750) * Math.sin(v_phi) * Math.sin(v_rho);
						v_z = (v_theta / (2 * pi)) * (750) * Math.cos(v_phi);
					} else if (which_one == 3) { //find phi as a function of theta and rho
						//find phi (rho = phi, theta = phi)
						v_rho = what_rho_is(v_phi, v_theta, which_one, scale);
						//spherical to rectangular
						v_x = (v_phi / pi) * (750) * Math.sin(v_rho) * Math.cos(v_theta);
						v_y = (v_phi / pi) * (750) * Math.sin(v_rho) * Math.sin(v_theta);
						v_z = (v_phi / pi) * (750) * Math.cos(v_rho);
					}
				} else if (what_system == 2) {//recgular coordinates
					//z = f(x, y)
					if (which_one == 1) {
						//change parameters
						//0 -> 2pi to -1000 -> 1000
						v_x = (v_theta * 750 / pi - 750);
						//0 -> pi to -1000 -> 1000
						v_y = (v_phi * 1500 / pi - 750);
						//find z
						v_z = scale * what_z_is(v_x / scale, v_y / scale, which_one);
					} else if (which_one == 2) { //y = f(x,z)
						//change parameters
						//0 -> 2pi to -1000 -> 1000
						v_x = (v_theta * 1000 / pi - 1000);
						//0 -> pi to -1000 -> 1000
						v_z = (v_phi * 2000 / pi - 1000);
						//find z
						v_y = scale * what_z_is(v_x / scale, v_z / scale, which_one);
					} else if (which_one == 3) { //x = f(z,y)
						//change parameters
						//0 -> 2pi to -1000 -> 1000
						v_z = (v_theta * 1000 / pi - 1000);
						//0 -> pi to -1000 -> 1000
						v_y = (v_phi * 2000 / pi - 1000);
						//find z
						v_x = scale * what_z_is(v_z / scale, v_y / scale, which_one);
					}
				} else if (what_system == 3) { //parametric
					//all 3 at once with pointers
					ParametricResult parametric_xyz = parametric_xyz(v_phi, v_theta);
					//blow up
					v_x = scale * parametric_xyz.x_uv;
					v_y = scale * parametric_xyz.y_uv;
					v_z = scale * parametric_xyz.z_uv;
				} else if (what_system == 4) { //cylindrical
					if (which_one == 1) {
						v_r = (v_phi * 500) / pi;
						//find z
						v_z = scale * what_polar_z_is(v_r, v_theta, which_one, scale) / scale;
						//find x, y
						v_x = scale * (v_r * Math.cos(v_theta));
						v_y = scale * (v_r * Math.sin(v_theta));
					} else if (which_one == 2) {
						//phi = r 0 -> pi to 0 -> 750
						v_z = (v_phi / pi) * (750) - (750 / 2);
						//find r
						v_r = scale * what_polar_z_is(v_z, v_theta, which_one, scale) / scale;
						//find x, y
						v_x = scale * (v_r * Math.cos(v_theta));
						v_y = scale * (v_r * Math.sin(v_theta));
					} else if (which_one == 3) {
						v_z = (v_theta / (2 * pi)) * (750) - (750 / 2);
						v_rho = (v_phi * 500) / pi;
						//find theta v_r = theta
						v_r = scale * what_polar_z_is(v_rho, v_z, which_one, scale) / scale;
						//find x, y
						v_x = scale * (v_rho) * Math.cos(v_r);
						v_y = scale * (v_rho) * Math.sin(v_r);
					}
				}

				if (v_z < 1920 && v_x < 1920 && v_y < 1920) {
					//find the 2d projection
					x_coord = determine_x_coord(v_x, v_y, v_z, t_rotation, p_rotation, s_rotation);
					y_coord = determine_y_coord(v_x, v_y, v_z, t_rotation, p_rotation, s_rotation);

					//find the disMath.tance to the plane
					c_dist = ((c1 * v_x) + (c2 * v_y) + (c3 * v_z) + c4) / c_norm;

					if ((x_coord >= 0) && (x_coord <= 1920) && (y_coord >= 0) && (y_coord <= 1200) && (is_it_closer(x_coord, y_coord, c_dist)) == true) {
						buffer.setColor(e_color(1 - Math.sin(v_phi)));
						buffer.drawLine(x_coord, y_coord, x_coord, y_coord);
					}
				} //end if (v_z)
			} //end for (theta)
			l_perc = (int) ((100 * v_phi) / pi);
			textprintf_ex(screen, 20, 20, makecol(255, 255, 255), String.format("loading %-5d...", l_perc));
		} //end for (phi)
	}

	//used to draw buttons for coordinate systems
	void draw_coord_buttons(int what_sys, int o_sys) {
		//redraw circles for what system sph rect para
		circle(buffer, 1500, 720, 15, makecol(255, 255, 255));
		circle(buffer, 1680, 720, 15, makecol(255, 255, 255));
		circle(buffer, 1850, 720, 15, makecol(255, 255, 255));
		circle(buffer, 1590, 720, 15, makecol(255, 255, 255));
		//circle(buffer, 1770, 720, 15, makecol(255,255,255));

		if (o_sys == 1) {
			circlefill(buffer, 1500, 720, 10, Color.BLACK);
		} else if (o_sys == 2) {
			circlefill(buffer, 1680, 720, 10, Color.BLACK);
		} else if (o_sys == 3) {
			circlefill(buffer, 1850, 720, 10, Color.BLACK);
		} else if (o_sys == 4) {
			circlefill(buffer, 1590, 720, 10, Color.BLACK);
			// if (o_sys == 5)
			//    circlefill(buffer, 1770, 720, 10, bitmap_mask_color(screen));
		}

		if (what_sys == 1) {
			circlefill(buffer, 1500, 720, 10, makecol(255, 255, 255));
		} else if (what_sys == 2) {
			circlefill(buffer, 1680, 720, 10, makecol(255, 255, 255));
		} else if (what_sys == 3) {
			circlefill(buffer, 1850, 720, 10, makecol(255, 255, 255));
		} else if (what_sys == 4) {
			circlefill(buffer, 1590, 720, 10, makecol(255, 255, 255));
			// if (what_sys == 5)
			//    circlefill(buffer, 1770, 720, 10,  makecol(255,255,255));
		}
	}

	///draws the interface based on the positions of the buttons
	void draw_interface(int mx, int mx2, int mx3, int mx4, int mx5, int w_sys, int o_w_sys, boolean dy) {
		//redraw bars detail/zoom/phi/theta/screen/redraw/rotate
		buffer.drawImage(bar, 1480, 780, null);
		buffer.drawImage(bar, 1480, 840, null);
		buffer.drawImage(bar, 1480, 900, null);
		buffer.drawImage(bar, 1480, 960, null);
		buffer.drawImage(bar, 1480, 1020, null);
		buffer.drawImage(rota2, 1730, 1100, null);

		//change button color
		if (dy) {
			buffer.drawImage(rdra2, 1530, 1100, null);
			dy = false;
		} else {
			buffer.drawImage(rdra, 1530, 1100, null);
			dy = true;
		}

		draw_coord_buttons(w_sys, o_w_sys);

		//redraw boxes on bars
		buffer.drawImage(box, mx4 - 8, 784, null);
		buffer.drawImage(box, mx5 - 8, 844, null);
		buffer.drawImage(box, mx3 - 8, 904, null);
		buffer.drawImage(box, mx - 8, 964, null);
		buffer.drawImage(box, mx2 - 8, 1040 - 16, null);
	}

	//w_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash)
	void draw_everything(int wh_s, int owh_s, double t_ro, double p_ro, double sc_ro, double de, double zo, int tr_b, int pr_b, int sr_b, int d_b, int s_b, boolean b_fl, int w_t) {
		//draw everything
		rectfill(buffer, 0, 0, bufferImage.getWidth(), bufferImage.getHeight(), Color.BLACK);
		draw_graph(wh_s, t_ro, p_ro, sc_ro, de, zo, w_t);
		draw_interface(tr_b, sr_b, pr_b, d_b, s_b, wh_s, owh_s, b_fl);
		draw_axes(t_ro, p_ro, sc_ro);
		textprintf_ex(buffer, 1490, 675, makecol(255, 255, 255), "sph        pol        rect                  par  ");
		textprintf_ex(buffer, 1490, 765, makecol(255, 255, 255), String.format("less          detail =  %-5f            more", de));
		textprintf_ex(buffer, 1490, 825, makecol(255, 255, 255), String.format("less        zoom    =    %-5f          more", 5 * Math.pow(10, (((double) (s_b - 1500) / 90) - 2))));
		textprintf_ex(buffer, 1490, 885, makecol(255, 255, 255), String.format("0             phi rotation =  %-5d           pi", (pr_b - 1500) / 2));
		textprintf_ex(buffer, 1490, 945, makecol(255, 255, 255), String.format("0            theta rotation = %-5d          2pi", tr_b - 1500));
		textprintf_ex(buffer, 1490, 1005, makecol(255, 255, 255), String.format("0           screen rotation = %-5d          2pi", sr_b - 1500));
		screen.drawImage(bufferImage, 0, 0, null);
	}

	private static class Graph3dMouseAdapter extends MouseAdapter {
		int mouse_x;
		int mouse_y;
		private final Graph3d graph;

		//rotations size and scale
		double detail = pi / 90;
		double zoom = 400;
		double th_rotation = pi / 4, ph_rotation = pi / 4, sc_rotation = 3 * pi / 2;
		//for button color switch
		boolean button_flash = true;
		//for displaying rotations ect
		int t_rot_bar_x = 1545, s_rot_bar_x = 1770, p_rot_bar_x = 1590, d_bar_x = 1500, s_bar_x = 1700;
		//what system 1 sph 2 rect 3 para and proir used system
		int w_s = 1, ol_s = 1, which_t = 1;

		public Graph3dMouseAdapter(Graph3d graph) {
			this.graph = graph;
		}

		@Override
		public void mouseMoved(MouseEvent e) {
			mouse_x = e.getX();
			mouse_y = e.getY();
		}

		@Override
		public void mouseDragged(MouseEvent e) {
			mouse_x = e.getX();
			mouse_y = e.getY();
			mouseDepressed();
		}

		@Override
		public void mousePressed(MouseEvent e) {
			//left click
			if (e.getButton() == MouseEvent.BUTTON1) {
				mouseDepressed();
			} else if (e.getButton() == MouseEvent.BUTTON3) { //right click
				which_t++;
				if (which_t >= 4) {
					which_t = 1;
				}
				graph.draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
			}
		}

		private void mouseDepressed() {
			//if the detail bar is clicked on
			if ((mouse_x >= 1500) && (mouse_x <= 1860) && (mouse_y >= 760) && (mouse_y <= 820)) {
				detail = mouse_x;
				d_bar_x = mouse_x;
				buffer.drawImage(bar, 1480, 780, null);
				buffer.drawImage(box, d_bar_x - 8, 784, null);
				//0 to pi/720
				detail = ((-7 * pi / 720) * ((detail - 1500) / 360)) + (pi / 90);

				//for fast rotation
				if (detail > 0.02) {
					graph.draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
				}
			} else if ((mouse_x >= 1500) && (mouse_x <= 1860) && (mouse_y >= 840) && (mouse_y <= 880)) { //if the size bar is clicked on
				zoom = mouse_x;
				s_bar_x = mouse_x;
				buffer.drawImage(bar, 1480, 840, null);
				buffer.drawImage(box, s_bar_x - 8, 844, null);

				//0 to 720
				zoom = 5 * Math.pow(10, ((double) (s_bar_x - 1500) / 90) - 2);

				//for fast rotation
				if (detail > 0.02) {
					graph.draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
				}
			} else if ((mouse_x >= 1500) && (mouse_x <= 1860) && (mouse_y >= 900) && (mouse_y <= 940)) { //if the phi rot bar is clicked on
				ph_rotation = mouse_x;
				p_rot_bar_x = mouse_x;

				buffer.drawImage(bar, 1480, 900, null);
				buffer.drawImage(box, p_rot_bar_x - 8, 904, null);

				//0 to pi
				ph_rotation = ((ph_rotation - 1500) * pi) / 360;

				//for fast rotation
				if (detail > 0.02) {
					graph.draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
				}
			} else if ((mouse_x >= 1500) && (mouse_x <= 1860) && (mouse_y >= 960) && (mouse_y <= 1000)) { //if the theta rot bar is clicked on
				th_rotation = mouse_x;
				t_rot_bar_x = mouse_x;

				buffer.drawImage(bar, 1480, 960, null);
				buffer.drawImage(box, t_rot_bar_x - 8, 964, null);

				//0 to 2pi
				th_rotation = ((th_rotation - 1500) * pi) / 180;

				//for fast rotation
				if (detail > 0.02) {
					graph.draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
				}
			} else if ((mouse_x >= 1500) && (mouse_x <= 1860) && (mouse_y >= 1020) && (mouse_y <= 1060)) { // if screen rotation is clicked
				sc_rotation = mouse_x;
				s_rot_bar_x = mouse_x;

				buffer.drawImage(bar, 1480, 1020, null);
				buffer.drawImage(box, s_rot_bar_x - 8, 1040 - 16, null);

				//0 to 2pi
				sc_rotation = ((sc_rotation - 1500) * pi) / 180;

				//for fast rotation
				if (detail > 0.02) {
					graph.draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
				}
			} else if ((mouse_x >= 1530) && (mouse_x <= 1630) && (mouse_y >= 1110) && (mouse_y <= 1160)) { //redraw
				//flash button
				if (button_flash) {
					button_flash = false;
				} else {
					button_flash = true;
				}
				//draw graph interface etc
				graph.draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
			} else if ((mouse_x >= 1730) && (mouse_x <= 1830) && (mouse_y >= 1110) && (mouse_y <= 1160)) { //go / rotate
				try {
					ImageIO.write(bufferImage, "bmp", new File(System.getProperty("user.home") + "\\Desktop\\graph.bmp"));
				} catch (IOException e1) {
					e1.printStackTrace();
				}
				if (button_flash) {
					buffer.drawImage(rota, 1730, 1100, null);
					button_flash = false;
				} else {
					buffer.drawImage(rota2, 1730, 1100, null);
					button_flash = true;
				}
			} else if (Math.pow(mouse_x - 1500, 2) + Math.pow(mouse_y - 720, 2) <= 225) { //spherical
				ol_s = w_s;
				w_s = 1;
				if (detail > 0.02) {
					graph.draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
				} else {
					graph.draw_coord_buttons(w_s, ol_s);
				}
			}
			if (Math.pow(mouse_x - 1680, 2) + Math.pow(mouse_y - 720, 2) <= 225) { //rectangular
				ol_s = w_s;
				w_s = 2;
				if (detail > 0.02) {
					graph.draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
				} else {
					graph.draw_coord_buttons(w_s, ol_s);
				}
			} else if (Math.pow(mouse_x - 1850, 2) + Math.pow(mouse_y - 720, 2) <= 225) { //parametric
				ol_s = w_s;
				w_s = 3;
				if (detail > 0.02) {
					graph.draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
				} else {
					graph.draw_coord_buttons(w_s, ol_s);
				}
			} else if (Math.pow(mouse_x - 1590, 2) + Math.pow(mouse_y - 720, 2) <= 225) { //polar
				ol_s = w_s;
				w_s = 4;
				if (detail > 0.02) {
					graph.draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
				} else {
					graph.draw_coord_buttons(w_s, ol_s);
				}
			} else if (Math.pow(mouse_x - 1770, 2) + Math.pow(mouse_y - 720, 2) <= 225) {//diff eq
				ol_s = w_s;
				w_s = 5;
				if (detail > 0.02) {
					graph.draw_everything(w_s, ol_s, th_rotation, ph_rotation, sc_rotation, detail, zoom, t_rot_bar_x, p_rot_bar_x, s_rot_bar_x, d_bar_x, s_bar_x, button_flash, which_t);
				} else {
					graph.draw_coord_buttons(w_s, ol_s);
				}
			}
		}
	}

	private static class Graph3dKeyAdapter extends KeyAdapter {

		@Override
		public void keyPressed(KeyEvent e) {
			if (e.getKeyCode() == KeyEvent.VK_ESCAPE) {
				System.exit(0);
			}
		}
	}

	public static void main(String[] args) {
		Graph3d graph = new Graph3d();

		JPanel mainPanel = new JPanel();

		mainPanel.setPreferredSize(new Dimension(SCREEN_W, SCREEN_H));
		Graph3dMouseAdapter ma = new Graph3dMouseAdapter(graph);
		mainPanel.addMouseListener(ma);
		mainPanel.addMouseMotionListener(ma);
		Graph3dKeyAdapter keyAdapter = new Graph3dKeyAdapter();
		mainPanel.addKeyListener(keyAdapter);

		JFrame mainFrame = new JFrame("3dGraph");
		mainFrame.setContentPane(mainPanel);
		mainFrame.pack();
		mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		mainFrame.setFocusable(false);
		mainFrame.setVisible(true);
		mainPanel.requestFocus();

		//palette
		//set_palette(desktop_palette);

		screen = mainPanel.getGraphics();
		//create a bitmap for double buffering equal to the screen size declared globally
		bufferImage = new BufferedImage(SCREEN_W, SCREEN_H, BufferedImage.TYPE_INT_RGB);
		buffer = bufferImage.createGraphics();

		//bitmap used for mouse cursor
		mouse_cursorImage = new BufferedImage(27, 27, BufferedImage.TYPE_INT_RGB);
		Graphics mouse_cursor = mouse_cursorImage.createGraphics();
		for (double c = 0; c < 10; c += .005) {
			circle(mouse_cursor, 13, 13, (int) c / 2, makecol(255 - (int) (255 - 16 * c), 0, 0));
		}

		//bar
		bar = new BufferedImage(400, 40, BufferedImage.TYPE_INT_RGB);
		Graphics barGraphics = bar.createGraphics();
		for (int k2 = 20; k2 < 200; k2++) {
			int r = (int) ((Math.sin((k2 - 20) * 3.15159 / 180) + 1) * 63);
			int g = (int) ((Math.sin((k2 - 20) * 3.15159 / 180) + 1) * 255);
			int b = (int) ((Math.sin((k2 - 20) * 3.15159 / 1800) + 1) * 255);
			circle(barGraphics, 400 - k2, 20, 5, makecol(r, g, b));
			circle(barGraphics, k2, 20, 5, makecol(r, g, b));
		}
		barGraphics.setColor(makecol(0, 200, 200));
		barGraphics.drawLine(20, 20, 380, 20);

		//nob
		box = new BufferedImage(16, 32, BufferedImage.TYPE_INT_RGB);
		Graphics2D boxGraphics = box.createGraphics();
		rectfill(boxGraphics, 0, 0, 16, 32, makecol(255, 63, 127));
		rectfill(boxGraphics, 3, 3, 12, 28, makecol(55, 255, 255));

		//redraw 1
		rdra = new BufferedImage(100, 60, BufferedImage.TYPE_INT_RGB);
		Graphics rdraGraphics = rdra.createGraphics();
		//redraw 2
		rdra2 = new BufferedImage(100, 60, BufferedImage.TYPE_INT_RGB);
		Graphics rdra2Graphics = rdra2.createGraphics();
		//rotate 1
		rota = new BufferedImage(100, 60, BufferedImage.TYPE_INT_RGB);
		Graphics rotaGraphics = rota.createGraphics();
		//rotate 2
		rota2 = new BufferedImage(100, 60, BufferedImage.TYPE_INT_RGB);
		Graphics rota2Graphics = rota2.createGraphics();

		//draw 2x2 buttons
		for (int k3 = 0; k3 < 20; k3++) {
			rect(rdraGraphics, k3, k3, 100 - k3, 60 - k3, makecol(55 - k3, k3 * 5, k3 + 75));
			textprintf_ex(rdraGraphics, 27, 27, makecol(123, 123, 231), "REDRAW");

			rect(rdra2Graphics, k3, k3, 100 - k3, 60 - k3, makecol(85 - k3, k3, 5 * k3));
			textprintf_ex(rdra2Graphics, 27, 27, makecol(123, 231, 123), "REDRAW");

			rect(rotaGraphics, k3, k3, 100 - k3, 60 - k3, makecol(55 - k3, k3 * 5, k3 + 75));
			textprintf_ex(rotaGraphics, 27, 27, makecol(123, 231, 123), " SAVE ");

			rect(rota2Graphics, k3, k3, 100 - k3, 60 - k3, makecol(55 - k3, k3 * 5, k3 + 75));
			textprintf_ex(rota2Graphics, 27, 27, makecol(123, 231, 123), " SAVE ");
		}

		graph.draw_everything(ma.w_s, ma.ol_s, ma.th_rotation, ma.ph_rotation, ma.sc_rotation, ma.detail, ma.zoom, ma.t_rot_bar_x, ma.p_rot_bar_x, ma.s_rot_bar_x, ma.d_bar_x, ma.s_bar_x,
				ma.button_flash, ma.which_t);
	}

	private static Color makecol(int r, int g, int b) {
		return new Color(Math.max(0, r) % 256, Math.max(0, g) % 256, Math.max(0, b) % 256);
	}

	private static void rect(Graphics g, int x0, int y0, int x1, int y1, Color c) {
		g.setColor(c);
		g.drawRect(x0, y0, x1 - x0, y1 - y0);
	}

	private static void rectfill(Graphics g, int x0, int y0, int x1, int y1, Color c) {
		g.setColor(c);
		g.fillRect(x0, y0, x1 - x0, y1 - y0);
	}

	private static void circle(Graphics g, int x, int y, int r, Color c) {
		g.setColor(c);
		g.drawOval(x - r, y - r, r * 2, r * 2);
	}

	private static void circlefill(Graphics g, int x, int y, int r, Color c) {
		g.setColor(c);
		g.fillOval(x - r, y - r, r * 2, r * 2);
	}

	private static void textprintf_ex(Graphics g, int x, int y, Color c, String str) {
		g.setFont(font);
		g.setColor(c);
		g.drawString(str, x, y + g.getFontMetrics().getHeight() / 2);
	}
}