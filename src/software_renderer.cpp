#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CS248 {


// Implements SoftwareRenderer //

// fill a sample location with color
void SoftwareRendererImp::fill_sample(int sx, int sy, const Color &color) {
  // Task 2: implement this function

  // check bounds of sample coords with width and height of ac
	if (sx < 0 || sx >= s_width) return;
	if (sy < 0 || sy >= s_height) return;

  Color pixel_color;
	float inv255 = 1.0 / 255.0;
	pixel_color.r = sample_buffer[4 * (sx + sy * s_width)] * inv255;
	pixel_color.g = sample_buffer[4 * (sx + sy * s_width) + 1] * inv255;
	pixel_color.b = sample_buffer[4 * (sx + sy * s_width) + 2] * inv255;
	pixel_color.a = sample_buffer[4 * (sx + sy * s_width) + 3] * inv255;

  pixel_color = ref->alpha_blending_helper(pixel_color, color);

  sample_buffer[4 * (sx + sy * s_width)] = (uint8_t)(pixel_color.r * 255);
  sample_buffer[4 * (sx + sy * s_width) + 1] = (uint8_t)(pixel_color.g * 255);
  sample_buffer[4 * (sx + sy * s_width) + 2] = (uint8_t)(pixel_color.b * 255);
  sample_buffer[4 * (sx + sy * s_width) + 3]= (uint8_t)(pixel_color.a * 255);
}

// fill samples in the entire pixel specified by pixel coordinates
void SoftwareRendererImp::fill_pixel(int x, int y, const Color &color) {

	// Task 2: Re-implement this function

	// check bounds
	// if (x < 0 || x >= width) return;
	// if (y < 0 || y >= height) return;

	// Color pixel_color;
	// float inv255 = 1.0 / 255.0;
	// pixel_color.r = pixel_buffer[4 * (x + y * width)] * inv255;
	// pixel_color.g = pixel_buffer[4 * (x + y * width) + 1] * inv255;
	// pixel_color.b = pixel_buffer[4 * (x + y * width) + 2] * inv255;
	// pixel_color.a = pixel_buffer[4 * (x + y * width) + 3] * inv255;

	// pixel_color = ref->alpha_blending_helper(pixel_color, color);

	// pixel_buffer[4 * (x + y * width)] = (uint8_t)(pixel_color.r * 255);
	// pixel_buffer[4 * (x + y * width) + 1] = (uint8_t)(pixel_color.g * 255);
	// pixel_buffer[4 * (x + y * width) + 2] = (uint8_t)(pixel_color.b * 255);
	// pixel_buffer[4 * (x + y * width) + 3] = (uint8_t)(pixel_color.a * 255);

  for (int i = 0; i < sqrt(sample_rate); i++) {
    for (int j = 0; j < sqrt(sample_rate); j++) {
      fill_sample(x * sqrt(sample_rate) + i, y * sqrt(sample_rate) + j, color);
    }
  }

}

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  transformation = canvas_to_screen;

  // canvas outline
  Vector2D a = transform(Vector2D(0, 0)); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width, 0)); b.x++; b.y--;
  Vector2D c = transform(Vector2D(0, svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width, svg.height)); d.x++; d.y++;

  svg_bbox_top_left = Vector2D(a.x+1, a.y+1);
  svg_bbox_bottom_right = Vector2D(d.x-1, d.y-1);

  // draw all elements
  for (size_t i = 0; i < svg.elements.size(); ++i) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to pixel buffer
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 2: 
  // You may want to modify this for supersampling support
  // freeing the buffer 
  this->sample_rate = sample_rate;
  this->s_width = this->width * sqrt(sample_rate);
  this->s_height = this->height * sqrt(sample_rate);

  if (this->sample_buffer) {
    free(this->sample_buffer);
  }
  // allocating for the sample buffer 
  this->sample_buffer = (unsigned char *)malloc(this->width * this->height * 4 * sample_rate);
  memset(this->sample_buffer, 255, 4 * this->s_width * this->s_height);
}

void SoftwareRendererImp::set_pixel_buffer( unsigned char* pixel_buffer,
                                             size_t width, size_t height ) {

  // Task 2: 
  // You may want to modify this for supersampling support
  this->pixel_buffer = pixel_buffer;
  this->width = width;
  this->height = height;

  this->s_width = width * sqrt(this->sample_rate);
  this->s_height = height * sqrt(this->sample_rate);
  if (this->sample_buffer) {
    free(this->sample_buffer);
  }
  this->sample_buffer = (unsigned char *)malloc(this->s_width * this->s_height * 4);
  memset(this->sample_buffer, 255, 4 * this->s_width * this->s_height);
}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

	// Task 3 (part 1):
	// Modify this to implement the transformation stack

	switch (element->type) {
	case POINT:
		draw_point(static_cast<Point&>(*element));
		break;
	case LINE:
		draw_line(static_cast<Line&>(*element));
		break;
	case POLYLINE:
		draw_polyline(static_cast<Polyline&>(*element));
		break;
	case RECT:
		draw_rect(static_cast<Rect&>(*element));
		break;
	case POLYGON:
		draw_polygon(static_cast<Polygon&>(*element));
		break;
	case ELLIPSE:
		draw_ellipse(static_cast<Ellipse&>(*element));
		break;
	case IMAGE:
		draw_image(static_cast<Image&>(*element));
		break;
	case GROUP:
		draw_group(static_cast<Group&>(*element));
		break;
	default:
		break;
	}

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Advanced Task
  // Implement ellipse rasterization

}

void SoftwareRendererImp::draw_image( Image& image ) {

  // Advanced Task
  // Render image element with rotation

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= width) return;
  if (sy < 0 || sy >= height) return;

  // fill sample - NOT doing alpha blending!
  // TODO: Call fill_pixel here to run alpha blending
  pixel_buffer[4 * (sx + sy * width)] = (uint8_t)(color.r * 255);
  pixel_buffer[4 * (sx + sy * width) + 1] = (uint8_t)(color.g * 255);
  pixel_buffer[4 * (sx + sy * width) + 2] = (uint8_t)(color.b * 255);
  pixel_buffer[4 * (sx + sy * width) + 3] = (uint8_t)(color.a * 255);

}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Task 0: 
  // Implement Bresenham's algorithm (delete the line below and implement your own)
  //ref->rasterize_line_helper(x0, y0, x1, y1, width, height, color, this);

  float dx = x1 - x0;
  if (dx < 0) {
    float temp = x1;
    x1 = x0;
    x0 = temp;

    temp = y1;
    y1 = y0;
    y0 = temp;

    dx *= -1;
  }

  float dy = y1 - y0;
  int eps = 0;
  float m = dy / dx;
  float eps_fl = 0.0;

  if (m >= -1 && m <= 1) {
    float y = y0;
    for (float x = x0; x <= x1; x++) {
      rasterize_point(x, y, color);

      // CASE: Positive Slope
      if (m >= 0.0) {
        eps += dy;
        // Multiplication by 2 can be implemented by left-shift 
        if ((eps << 1) >= dx) {
          y++;  
          eps -= dx;
        }

      // CASE: Negative Slope
      } else {
        eps_fl += m;
        if (eps_fl <= -0.5) {
          y--;
          eps_fl++;
        }
      }
    }
  } else {
    if (dy < 0) {
      float temp = x1;
      x1 = x0;
      x0 = temp;

      temp = y1;
      y1 = y0;
      y0 = temp;

      dy *= -1;
      dx *= -1;
      m = dx / dy;
    }
    float x = x0;

    for (float y = y0; y <= y1; y++) {
      rasterize_point(x, y, color);

      // CASE: Positive Slope
      if (m >= 0.0) {
        eps += dx;
        // Multiplication by 2 can be implemented by left-shift 
        if ((eps << 1) >= dy) {
          x++;  
          eps -= dy;
        }

      // CASE: Negative Slope
      } else {
        eps_fl += m;
        if (eps_fl <= -0.5) {
          x--;
          eps_fl++;
        }
      }
    }
  }

  // Advanced Task
  // Drawing Smooth Lines with Line Width
}

std::vector<float> SoftwareRendererImp::compute_line_coefficients(float x0, float y0, float x1, float y1) {
  std::vector<float> line_coefficients;
  float A = y1 - y0;
  float B = x0 - x1;
  float C = (y0 * (x1 - x0)) - (x0 * (y1 - y0));
  line_coefficients.push_back(A);
  line_coefficients.push_back(B);
  line_coefficients.push_back(C);

  return line_coefficients;
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // compute boundary box
  float xmin = x0 < x1? (x0 < x2 ? x0 : x2) : (x1 < x2 ? x1 : x2);
  float xmax = x0 > x1? (x0 > x2 ? x0 : x2) : (x1 > x2 ? x1 : x2);
  float ymin = y0 < y1? (y0 < y2 ? y0 : y2) : (y1 < y2 ? y1 : y2);
  float ymax = y0 > y1? (y0 > y2 ? y0 : y2) : (y1 > y2 ? y1 : y2);

  // if vertices are clock-wise, switch points 1 and 2
  float orientation = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
  if (orientation < 0) {
    float temp = x1;
    x1 = x2;
    x2 = temp;
    temp = y1;
    y1 = y2;
    y2 = temp;
  }

  // compute coefficients for each line 
  std::vector<float> l1 = compute_line_coefficients(x0, y0, x1, y1);
  std::vector<float> l2 = compute_line_coefficients(x1, y1, x2, y2);
  std::vector<float> l3 = compute_line_coefficients(x2, y2, x0, y0);

  // for all pixels in bounding box
  for (int sx = (int)floor(xmin); sx <= (int)floor(xmax); sx++) {
    for (int sy = (int)floor(ymin); sy <= (int)floor(ymax); sy++) {
      // compute pixel center
      float xcenter = (float)sx + 0.5f;
      float ycenter = (float)sy + 0.5f;

      // determine whether pixel center is in triangle by checking each line
      if (((l1[0] * xcenter) + (l1[1] * ycenter) + l1[2]) <= 0) {
        if (((l2[0] * xcenter) + (l2[1] * ycenter) + l2[2]) <= 0) {
          if (((l3[0] * xcenter) + (l3[1] * ycenter) + l3[2]) <= 0) {
            fill_pixel(sx, sy, color);
          }
        }
      }
    }
  }

  // Advanced Task
  // Implementing Triangle Edge Rules

}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 4: 
  // Implement image rasterization

}

// resolve samples to pixel buffer
void SoftwareRendererImp::resolve( void ) {

  // Task 2: 
  // Implement supersampling

  // iterate over every pixel and "average" all the samples
  for (int x = 0; x < width; x++) {     // x coord
    for (int y = 0; y < height; y++) {  // y coord
      for (int k = 0; k < 4; k++) {     // color offset
        float color = 0;
        for (int i = 0; i < sqrt(sample_rate); i++) {    // x sample offset
          for (int j = 0; j < sqrt(sample_rate); j++) {  // y sample offset
            int sx = x * sqrt(sample_rate) + i;
            int sy = y * sqrt(sample_rate) + j;
            color += sample_buffer[4 * (sx + sy * s_width) + k];
          }
        }
        color /= sample_rate;
        pixel_buffer[4 * (x + y * width) + k] = color;
      }
    }
  }

}

Color SoftwareRendererImp::alpha_blending(Color pixel_color, Color color)
{
  // Task 5
  // Implement alpha compositing
  return pixel_color;
}


} // namespace CS248
