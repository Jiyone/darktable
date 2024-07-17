/*
    This file is part of darktable,
    Copyright (C) 2020 darktable developers.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "bauhaus/bauhaus.h"
#include "chart/common.h"
#include "dtgtk/drawingarea.h"
#include "common/darktable.h"
#include "common/opencl.h"
#include "common/colorchecker.h"
#include "chart/colorchart.h"
#include "control/control.h"
#include "develop/develop.h"
#include "develop/imageop.h"
#include "develop/imageop_math.h"
#include "develop/imageop_gui.h"
#include "develop/openmp_maths.h"
#include "dtgtk/button.h"
#include "dtgtk/resetlabel.h"
#include "gui/accelerators.h"
#include "gui/gtk.h"
#include "gui/presets.h"
#include "gui/color_picker_proxy.h"
#include "iop/iop_api.h"

#include <glib.h>
#include <math.h>
#include <stdlib.h>
#include <gtk/gtk.h>

#if defined(__GNUC__)
#pragma GCC optimize ("unroll-loops", "tree-loop-if-convert", \
                      "tree-loop-distribution", "no-strict-aliasing", \
                      "loop-interchange", "loop-nest-optimize", "tree-loop-im", \
                      "unswitch-loops", "tree-loop-ivcanon", "ira-loop-pressure", \
                      "split-ivs-in-unroller", "variable-expansion-in-unroller", \
                      "split-loops", "ivopts", "predictive-commoning",\
                      "tree-loop-linear", "loop-block", "loop-strip-mine", \
                      "finite-math-only", "fp-contract=fast", "fast-math", \
                      "tree-vectorize", "no-math-errno")
#endif

/** DOCUMENTATION
 *
 * This module allows to invert scanned negatives and simulate their print on paper, based on Kodak Cineon
 * densitometry algorithm. It is better than the old invert module because it takes into account the Dmax of the film
 * and allows white balance adjustments, as well as paper grade (gamma) simulation. It also allows density correction
 * in log space, to account for the exposure settings of the scanner. Finally, it is applied after input colour profiling,
 * which means the inversion happens after the scanner or the camera got color-corrected, while the old invert module
 * invert the RAW, non-demosaiced, file before any colour correction.
 *
 * References :
 *
 *  - https://www.kodak.com/uploadedfiles/motion/US_plugins_acrobat_en_motion_education_sensitometry_workbook.pdf
 *  - http://www.digital-intermediate.co.uk/film/pdf/Cineon.pdf
 *  - https://lists.gnu.org/archive/html/openexr-devel/2005-03/msg00009.html
 **/

 #define THRESHOLD 2.3283064365386963e-10f // -32 EV


DT_MODULE_INTROSPECTION(2, dt_iop_negadoctor_params_t)


typedef enum dt_iop_negadoctor_filmstock_t
{
  // What kind of emulsion are we working on ?
  DT_FILMSTOCK_NB = 0,   // $DESCRIPTION: "black and white film"
  DT_FILMSTOCK_COLOR = 1 // $DESCRIPTION: "color film"
} dt_iop_negadoctor_filmstock_t;

typedef enum density_val_t
{
  RANGE_MAX_RED = 0,
  RANGE_MAX_GREEN,
  RANGE_MAX_BLUE,

  RANGE_MIN_RED,
  RANGE_MIN_GREEN,
  RANGE_MIN_BLUE,

  SPOT_RED,
  SPOT_GREEN,
  SPOT_BLUE,

  RANGE_MAX_AVERAGE,
  RANGE_MIN_AVERAGE,
  SPOT_AVERAGE,

  RANGE_LAST_FIELD
} density_val_t;

typedef struct dt_extraction_data_t
{
  int width;
  int height;
  float radius_x;
  float radius_y;
  float *homography;
  float *inverse_homography;
}dt_extraction_data_t;

typedef struct dt_iop_negadoctor_params_t
{
  dt_iop_negadoctor_filmstock_t film_stock; /* $DEFAULT: DT_FILMSTOCK_COLOR $DESCRIPTION: "film stock" */
  float Dmin[4];                            /* color of film substrate
                                               $MIN: 0.00001 $MAX: 1.5 $DEFAULT: 1.0 */
  float wb_high[4];                         /* white balance RGB coeffs (illuminant)
                                               $MIN: 0.25 $MAX: 2 $DEFAULT: 1.0 */
  float wb_low[4];                          /* white balance RGB offsets (base light)
                                               $MIN: -2.0 $MAX: 2.0 $DEFAULT: 0.0 */
  float D_max;                              /* max density of film
                                               $MIN: 0.1 $MAX: 6 $DEFAULT: 2.046 */
  float offset;                             /* inversion offset
                                               $MIN: -1.0 $MAX: 1.0 $DEFAULT: -0.00 $DESCRIPTION: "scan exposure bias" */
  float D[RANGE_LAST_FIELD];                /* Densitometer reading
                                               $MIN: 0.01 $DEFAULT: 0.01 $DESCRIPTION: "Density measurments" */
  int checker;                              /* Step wedges
                                               $DEFAULT: 0 $DESCRIPTION: "Chosen density chart" */
  float D_film_max[3];                      /* Film's Dmax input
                                               $MIN: 0.0 $MAX: 4.0 $DEFAULT: 2.98 $DESCRIPTION: "Enter the film's dmax (from film's datasheet)" */
  float black;                              /* display black level
                                               $MIN: -0.5 $MAX: 0.5 $DEFAULT: 0.0 $DESCRIPTION: "paper black (density correction)" */
  float gamma;                              /* display gamma
                                               $MIN: 1.0 $MAX: 8.0 $DEFAULT: 4.0 $DESCRIPTION: "paper grade (gamma)" */
  float soft_clip;                          /* highlights roll-off
                                               $MIN: 0.0001 $MAX: 1.0 $DEFAULT: 1.0 $DESCRIPTION: "paper gloss (specular highlights)" */
  float exposure;                           /* extra exposure
                                               $MIN: -2.0 $MAX: 2.0 $DEFAULT: 1.0 $DESCRIPTION: "print exposure adjustment" */
} dt_iop_negadoctor_params_t;


typedef struct dt_iop_negadoctor_data_t
{
  dt_aligned_pixel_t Dmin;                // color of film substrate
  dt_aligned_pixel_t wb_high;             // white balance RGB coeffs / Dmax
  dt_aligned_pixel_t offset;              // inversion offset
  float black;                            // display black level
  float gamma;                            // display gamma
  float soft_clip;                        // highlights roll-off
  float soft_clip_comp;                   // 1 - softclip, complement to 1
  float exposure;                         // extra exposure
} dt_iop_negadoctor_data_t;


typedef struct dt_iop_negadoctor_gui_data_t
{
  GtkNotebook *notebook;
  GtkWidget *film_stock;
  GtkWidget *Dmin_R, *Dmin_G, *Dmin_B;
  GtkWidget *wb_high_R, *wb_high_G, *wb_high_B;
  GtkWidget *wb_low_R, *wb_low_G, *wb_low_B;
  GtkWidget *D_max;
  GtkWidget *offset;

  dt_gui_collapsible_section_t cs;

  GtkWidget *density_info[RANGE_LAST_FIELD];
  GtkWidget *D_sampler, *Dmin_range_sampler, *Dmax_range_sampler;
  GtkWidget *D_film_max_R, *D_film_max_G, *D_film_max_B;
  GtkWidget *checkers_list, *safety;
  GtkWidget *button_commit;
  gboolean run_validation;      // order a profile validation at next pipeline recompute
  float homography[9];          // the perspective correction matrix
  float inverse_homography[9];  // The inverse perspective correction matrix
  gboolean checker_ready;       // notify that a checker bounding box is ready to be used
  
  gboolean is_profiling_started;
  point_t box[4];           // the current coordinates, possibly non rectangle, of the bounding box for the color checker
  point_t ideal_box[4];     // the desired coordinates of the perfect rectangle bounding box for the color checker
  point_t center_box;       // the barycenter of both boxes
  gboolean active_node[4];  // true if the cursor is close to a node (node = corner of the bounding box)
  gboolean is_cursor_close; // do we have the cursor close to a node ?
  gboolean drag_drop;       // are we currently dragging and dropping a node ?
  point_t click_start;      // the coordinates where the drag and drop started
  point_t click_end;        // the coordinates where the drag and drop started
  dt_step_wedge_t *checker;
  float safety_margin;
  
  GtkWidget *black, *gamma, *soft_clip, *exposure;
  GtkWidget *Dmin_picker, *Dmin_sampler;
  GtkWidget *WB_high_picker, *WB_high_norm, *WB_high_sampler;
  GtkWidget *WB_low_picker, *WB_low_norm, *WB_low_sampler;
} dt_iop_negadoctor_gui_data_t;


typedef struct dt_iop_negadoctor_global_data_t
{
  int kernel_negadoctor;
} dt_iop_negadoctor_global_data_t;

typedef struct density_result_t
{
  float ref_Dmax;
  float ref_Dmin;
  float Dmax_test[3];
  float Dmin_test[3];
} density_result_t;

const char *name()
{
  return _("negadoctor");
}

const char *aliases()
{
  return _("film|invert|negative|scan");
}

const char **description(struct dt_iop_module_t *self)
{
  return dt_iop_set_description(self, _("invert film negative scans and simulate printing on paper"),
                                      _("corrective and creative"),
                                      _("linear, RGB, display-referred"),
                                      _("non-linear, RGB"),
                                      _("non-linear, RGB, display-referred"));
}

int flags()
{
  return IOP_FLAGS_INCLUDE_IN_STYLES | IOP_FLAGS_ALLOW_TILING | IOP_FLAGS_ONE_INSTANCE;
}


int default_group()
{
  return IOP_GROUP_FILM;
}


int default_colorspace(dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  return IOP_CS_RGB;
}

int legacy_params(dt_iop_module_t *self, const void *const old_params, const int old_version,
                  void *new_params, const int new_version)
{
  if(old_version == 1 && new_version == 2)
  {
    typedef struct dt_iop_negadoctor_params_v1_t
    {
      dt_iop_negadoctor_filmstock_t film_stock;
      dt_aligned_pixel_t Dmin;                // color of film substrate
      dt_aligned_pixel_t wb_high;             // white balance RGB coeffs (illuminant)
      dt_aligned_pixel_t wb_low;              // white balance RGB offsets (base light)
      float D_max;                            // max density of film
      float offset;                           // inversion offset
      float black;                            // display black level
      float gamma;                            // display gamma
      float soft_clip;                        // highlights roll-off
      float exposure;                         // extra exposure
    } dt_iop_negadoctor_params_v1_t;

    dt_iop_negadoctor_params_v1_t *o = (dt_iop_negadoctor_params_v1_t *)old_params;
    dt_iop_negadoctor_params_t *n = (dt_iop_negadoctor_params_t *)new_params;
    dt_iop_negadoctor_params_t *d = (dt_iop_negadoctor_params_t *)self->default_params;

    *n = *d; // start with a fresh copy of default parameters

    // WARNING: when copying the arrays in a for loop, gcc wrongly assumed
    //          that n and o were aligned and used AVX instructions for me,
    //          which segfaulted. let's hope this doesn't get optimized too much.
    n->film_stock = o->film_stock;
    n->Dmin[0] = o->Dmin[0];
    n->Dmin[1] = o->Dmin[1];
    n->Dmin[2] = o->Dmin[2];
    n->Dmin[3] = o->Dmin[3];
    n->wb_high[0] = o->wb_high[0];
    n->wb_high[1] = o->wb_high[1];
    n->wb_high[2] = o->wb_high[2];
    n->wb_high[3] = o->wb_high[3];
    n->wb_low[0] = o->wb_low[0];
    n->wb_low[1] = o->wb_low[1];
    n->wb_low[2] = o->wb_low[2];
    n->wb_low[3] = o->wb_low[3];
    n->D_max = o->D_max;
    n->offset = o->offset;
    n->black = o->black;
    n->gamma = o->gamma;
    n->soft_clip = o->soft_clip;
    n->exposure = o->exposure;

    return 0;
  }
  return 1;
}

static inline void _get_densities(float *patches, dt_aligned_pixel_t pix)
{
  for(int c = 0; c < 3; c++) pix[c] = log10f(1.f / patches[c]);
}

static void _get_patch_average(const float *const restrict in, const dt_step_wedge_t *checker, float *const restrict patches, const dt_extraction_data_t dim)
{
  /* Get the average color over each patch */
  for(size_t k = 0; k < checker->patches; k++)
  {
    // center of the patch in the ideal reference
    const point_t center = { checker->values[k].x, checker->values[k].y };

    // corners of the patch in the ideal reference
    const point_t corners[4] = { {center.x - dim.radius_x, center.y - dim.radius_y},
                                 {center.x + dim.radius_x, center.y - dim.radius_y},
                                 {center.x + dim.radius_x, center.y + dim.radius_y},
                                 {center.x - dim.radius_x, center.y + dim.radius_y} };

    // apply patch coordinates transform depending on perspective
    point_t new_corners[4];
    // find the bounding box of the patch at the same time
    size_t x_min = dim.width - 1;
    size_t x_max = 0;
    size_t y_min = dim.height - 1;
    size_t y_max = 0;
    for(size_t c = 0; c < 4; c++) {
      new_corners[c] = apply_homography(corners[c], dim.homography);
      x_min = fminf(new_corners[c].x, x_min);
      x_max = fmaxf(new_corners[c].x, x_max);
      y_min = fminf(new_corners[c].y, y_min);
      y_max = fmaxf(new_corners[c].y, y_max);
    }

    x_min = CLAMP((size_t)floorf(x_min), 0, dim.width - 1);
    x_max = CLAMP((size_t)ceilf(x_max), 0, dim.width - 1);
    y_min = CLAMP((size_t)floorf(y_min), 0, dim.height - 1);
    y_max = CLAMP((size_t)ceilf(y_max), 0, dim.height - 1);

    //// Get the average color on the patch
    patches[k * 4] = patches[k * 4 + 1] = patches[k * 4 + 2] = patches[k * 4 + 3] = 0.f;
    size_t num_elem = 0;

    // Loop through the rectangular bounding box
    for(size_t j = y_min; j < y_max; j++)
      for(size_t i = x_min; i < x_max; i++)
      {
        // Check if this pixel lies inside the sampling area and sample if it does
        point_t current_point = { i + 0.5f, j + 0.5f };
        current_point = apply_homography(current_point, dim.inverse_homography);
        current_point.x -= center.x;
        current_point.y -= center.y;

        if(current_point.x < dim.radius_x && current_point.x > -dim.radius_x &&
           current_point.y < dim.radius_y && current_point.y > -dim.radius_y)
        {
          for(size_t c = 0; c < 3; c++)
          {
            patches[k * 4 + c] += in[(j * dim.width + i) * 4 + c];

            // Debug : inpaint a black square in the preview to ensure the coordinates of
            // overlay drawings and actual pixel processing match
            // out[(j * width + i) * 4 + c] = 0.f;
          }
          num_elem++;
        }
      }

    for(size_t c = 0; c < 3; c++) patches[k * 4 + c] /= (float)num_elem;

    // Convert to density
    
    dt_aligned_pixel_t densities;
    _get_densities(patches + k * 4, densities);
    for(size_t o = 0; o < 3; o++) patches[k * 4 + o] = densities[o];

    fprintf(stdout, "[NEG ref] Patch %i: %.2f %.2f %.2f\n", k, patches[k * 4], patches[k * 4 + 1], patches[k * 4 + 2]);
  }
  fprintf(stdout, "[NEG ref] End of _get_patch_average\n\n");
}

static const density_result_t _extract_patches(const float *const restrict in, const dt_iop_roi_t *const roi_in, dt_iop_negadoctor_gui_data_t *g,
                                                  float *const restrict patches)
{
  const float rx = g->checker->radius * hypotf(1.f, g->checker->ratio) * g->safety_margin;

  dt_extraction_data_t dimensions = { roi_in->width,
                                            roi_in->height,
                                            rx,
                                            rx / g->checker->ratio,
                                            g->homography,
                                            g->inverse_homography };

  _get_patch_average(in, g->checker, patches, dimensions);

  // find reference white and black patch
  float ref_Dmax = g->checker->values[g->checker->dmax].density;
  float ref_Dmin = g->checker->values[g->checker->dmin].density;

  fprintf(stdout, "[NEG ref] Dmax: %f\t Dmin: %f\n", ref_Dmax, ref_Dmin);
  // find test dmax and dmin patch
  float Dmax_test[3];
  for(size_t c = 0; c < 3; c++)
  {
    Dmax_test[c] = patches[g->checker->dmax * 4 + c];
    fprintf(stdout, "[NEG ref] Dmax_test: %f\n", Dmax_test[c]);
  }

  float Dmin_test[3];
  for(size_t c = 0; c < 3; c++)
  {
    Dmin_test[c] = patches[g->checker->dmin * 4 + c];
    fprintf(stdout, "[NEG ref] Dmin_test: %f\n", Dmin_test[c]);
  }

  const density_result_t result = { ref_Dmax, ref_Dmin,
                                    { Dmax_test[0], Dmax_test[1], Dmax_test[2] },
                                    { Dmin_test[0], Dmin_test[1], Dmin_test[2] } };
  return result;
}

void validate_color_checker(const float *const restrict in, dt_dev_pixelpipe_iop_t *const piece,
                            const dt_iop_roi_t *const roi_in, struct dt_iop_module_t *const self)
{
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;
  //const dt_iop_negadoctor_data_t *const d = piece->data;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;

  float *const restrict patches = dt_alloc_sse_ps(4 * g->checker->patches);
  density_result_t val = _extract_patches(in, roi_in, g, patches);

  // Update GUI label

  for(int c = 0; c < 3; c++)
  {
    p->D[c] = val.Dmax_test[c];
    p->D[c + 3] = val.Dmax_test[c + 3];
    p->D_film_max[c] = val.ref_Dmax;
  }
  fprintf(stdout, "[NEG ref] range D max: %f %f %f\n", p->D[0], p->D[1], p->D[2]);
  fprintf(stdout, "[NEG ref] range D min: %f %f %f\n", p->D[3], p->D[4], p->D[5]);
  fprintf(stdout, "[NEG ref] ref dmax: %f\n", p->D[0], p->D_film_max);

  float thing[3] = { 0.f };
  float average[2][3] = {{ 0.f }};

  for(int k = 0; k < 9; k++) 
  {
    size_t channel = k % 3;

    thing[channel] = p->D[k] * val.ref_Dmax / val.Dmax_test[channel];
    //g_free(g->density_info[k]);
    gtk_label_set_text(GTK_LABEL(g->density_info[k]), g_strdup_printf("%.2f", thing[channel]));

    if(k % 3 == 2)
    {
      average[0][k / 3] = v_minf(thing);
      average[1][k / 3] = v_maxf(thing);
    }
  }

  for(int a = 0; a < 3; a++)
  {
    float averaged = (average[0][a] + average[1][a]) / 2;
    //g_free(g->density_info[a + 9]);
    gtk_label_set_text(GTK_LABEL(g->density_info[a + 9]), g_strdup_printf("%.2f", averaged));
  }

  dt_free_align(patches);
}

void commit_params(dt_iop_module_t *self, dt_iop_params_t *p1, dt_dev_pixelpipe_t *pipe,
                   dt_dev_pixelpipe_iop_t *piece)
{
  const dt_iop_negadoctor_params_t *const p = (dt_iop_negadoctor_params_t *)p1;
  dt_iop_negadoctor_data_t *const d = (dt_iop_negadoctor_data_t *)piece->data;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;

  // keep WB_high even in B&W mode to apply sepia or warm tone look
  // but premultiply it aheard with Dmax to spare one div per pixel
  for(size_t c = 0; c < 4; c++) d->wb_high[c] = p->wb_high[c] / p->D_max;

  for(size_t c = 0; c < 4; c++) d->offset[c] = p->wb_high[c] * p->offset + p->wb_low[c] * -0.3f;

  // ensure we use a monochrome Dmin for B&W film
  if(p->film_stock == DT_FILMSTOCK_COLOR)
    for(size_t c = 0; c < 4; c++) d->Dmin[c] = p->Dmin[c];
  else if(p->film_stock == DT_FILMSTOCK_NB)
    for(size_t c = 0; c < 4; c++) d->Dmin[c] = p->Dmin[0];

  // arithmetic trick allowing to rewrite the pixel inversion as FMA
  d->black = -p->exposure * (1.0f + p->black);

  // highlights soft clip
  d->soft_clip = p->soft_clip;
  d->soft_clip_comp = 1.0f - p->soft_clip;

  // copy
  d->exposure = p->exposure;
  d->gamma = p->gamma;

  // Disable OpenCL path if we are in any kind of diagnose mode (only C path has diagnostics)
  if(self->dev->gui_attached && g)
  {
    if( /*(g->run_profile && piece->pipe->type == DT_DEV_PIXELPIPE_PREVIEW) || // color checker extraction mode*/
        (g->run_validation && piece->pipe->type == DT_DEV_PIXELPIPE_PREVIEW) /*|| // delta E validation
        ( (d->illuminant_type == DT_ILLUMINANT_DETECT_EDGES ||
           d->illuminant_type == DT_ILLUMINANT_DETECT_SURFACES ) && // WB extraction mode
           piece->pipe->type == DT_DEV_PIXELPIPE_FULL )*/ )
    {
      piece->process_cl_ready = 0;
    }
  }
}

void process(struct dt_iop_module_t *const self, dt_dev_pixelpipe_iop_t *const piece,
             const void *const restrict ivoid, void *const restrict ovoid,
             const dt_iop_roi_t *const restrict roi_in, const dt_iop_roi_t *const restrict roi_out)
{
  const dt_iop_negadoctor_data_t *const d = piece->data;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  assert(piece->colors = 4);

  const float *const restrict in = (float *)ivoid;
  float *const restrict out = (float *)ovoid;


#ifdef _OPENMP
  #pragma omp parallel for simd default(none) \
    dt_omp_firstprivate(d, in, out, roi_out) \
    aligned(in, out:64) collapse(2)
#endif
  for(size_t k = 0; k < (size_t)roi_out->height * roi_out->width * 4; k += 4)
  {
    for(size_t c = 0; c < 4; c++)
    {
      // Unpack vectors one by one with extra pragmas to be sure the compiler understands they can be vectorized
      const float *const restrict pix_in = in + k;
      float *const restrict pix_out = out + k;
      const float *const restrict Dmin = __builtin_assume_aligned(d->Dmin, 16);
      const float *const restrict wb_high = __builtin_assume_aligned(d->wb_high, 16);
      const float *const restrict offset = __builtin_assume_aligned(d->offset, 16);

      // Convert transmission to density using Dmin as a fulcrum
      const float density = - log10f(Dmin[c] / fmaxf(pix_in[c], THRESHOLD)); // threshold to -32 EV

      // Correct density in log space
      const float corrected_de = wb_high[c] * density + offset[c];

      // Print density on paper : ((1 - 10^corrected_de + black) * exposure)^gamma rewritten for FMA
      const float print_linear = -(d->exposure * fast_exp10f(corrected_de) + d->black);
      const float print_gamma = powf(fmaxf(print_linear, 0.0f), d->gamma); // note : this is always > 0

      // Compress highlights. from https://lists.gnu.org/archive/html/openexr-devel/2005-03/msg00009.html
      pix_out[c] =  (print_gamma > d->soft_clip) ? d->soft_clip + (1.0f - fast_expf(-(print_gamma - d->soft_clip) / d->soft_clip_comp)) * d->soft_clip_comp
                                                 : print_gamma;
    }
  }

  if(piece->pipe->mask_display & DT_DEV_PIXELPIPE_DISPLAY_MASK)
    dt_iop_alpha_copy(ivoid, ovoid, roi_out->width, roi_out->height);

  // Density reading
  if(self->dev->gui_attached && g)
  {  
    fprintf(stdout, "[NEG ref] run?: %s\n",g->run_validation ? "YES" : "NO" );
    fprintf(stdout, "[NEG ref] checker: %i\n", g->checker);
    if(g->run_validation && piece->pipe->type == DT_DEV_PIXELPIPE_PREVIEW && g->checker)
    {
      validate_color_checker(in, piece, roi_out, self);
      g->run_validation = FALSE;
    }
  }
}


#ifdef HAVE_OPENCL
int process_cl(struct dt_iop_module_t *const self, dt_dev_pixelpipe_iop_t *const piece, cl_mem dev_in, cl_mem dev_out,
               const dt_iop_roi_t *const restrict roi_in, const dt_iop_roi_t *const restrict roi_out)
{
  const dt_iop_negadoctor_data_t *const d = (dt_iop_negadoctor_data_t *)piece->data;
  const dt_iop_negadoctor_global_data_t *const gd = (dt_iop_negadoctor_global_data_t *)self->global_data;

  cl_int err = -999;

  const int devid = piece->pipe->devid;
  const int width = roi_in->width;
  const int height = roi_in->height;

  size_t sizes[] = { ROUNDUPDWD(width, devid), ROUNDUPDHT(height, devid), 1 };

  dt_opencl_set_kernel_arg(devid, gd->kernel_negadoctor, 0, sizeof(cl_mem), (void *)&dev_in);
  dt_opencl_set_kernel_arg(devid, gd->kernel_negadoctor, 1, sizeof(cl_mem), (void *)&dev_out);
  dt_opencl_set_kernel_arg(devid, gd->kernel_negadoctor, 2, sizeof(int), (void *)&width);
  dt_opencl_set_kernel_arg(devid, gd->kernel_negadoctor, 3, sizeof(int), (void *)&height);
  dt_opencl_set_kernel_arg(devid, gd->kernel_negadoctor, 4, 4 * sizeof(float), (void *)&d->Dmin);
  dt_opencl_set_kernel_arg(devid, gd->kernel_negadoctor, 5, 4 * sizeof(float), (void *)&d->wb_high);
  dt_opencl_set_kernel_arg(devid, gd->kernel_negadoctor, 6, 4 * sizeof(float), (void *)&d->offset);
  dt_opencl_set_kernel_arg(devid, gd->kernel_negadoctor, 7, sizeof(float), (void *)&d->exposure);
  dt_opencl_set_kernel_arg(devid, gd->kernel_negadoctor, 8, sizeof(float), (void *)&d->black);
  dt_opencl_set_kernel_arg(devid, gd->kernel_negadoctor, 9, sizeof(float), (void *)&d->gamma);
  dt_opencl_set_kernel_arg(devid, gd->kernel_negadoctor, 10, sizeof(float), (void *)&d->soft_clip);
  dt_opencl_set_kernel_arg(devid, gd->kernel_negadoctor, 11, sizeof(float), (void *)&d->soft_clip_comp);

  err = dt_opencl_enqueue_kernel_2d(devid, gd->kernel_negadoctor, sizes);
  if(err != CL_SUCCESS) goto error;
    return TRUE;

error:
  dt_print(DT_DEBUG_OPENCL, "[opencl_negadoctor] couldn't enqueue kernel! %d\n", err);
  return FALSE;
}
#endif


void init(dt_iop_module_t *module)
{
  dt_iop_default_init(module);

  dt_iop_negadoctor_params_t *d = module->default_params;

  d->Dmin[0] = 1.00f;
  d->Dmin[1] = 0.45f;
  d->Dmin[2] = 0.25f;
}

void init_presets(dt_iop_module_so_t *self)
{
  dt_iop_negadoctor_params_t tmp = (dt_iop_negadoctor_params_t){ .film_stock = DT_FILMSTOCK_COLOR,
                                                                 .Dmin = { 1.13f, 0.49f, 0.27f, 0.0f},
                                                                 .wb_high = { 1.0f, 1.0f, 1.0f, 0.0f },
                                                                 .wb_low = { 1.0f, 1.0f, 1.0f, 0.0f },
                                                                 .D_max = 1.6f,
                                                                 .offset = -0.05f,
                                                                 .gamma = 4.0f,
                                                                 .soft_clip = 0.75f,
                                                                 .exposure = 0.9245f,
                                                                 .black = 0.0755f,
                                                                 .checker = 0};


  dt_gui_presets_add_generic(_("color film"), self->op,
                             self->version(), &tmp, sizeof(tmp), 1, DEVELOP_BLEND_CS_RGB_DISPLAY);

  dt_iop_negadoctor_params_t tmq = (dt_iop_negadoctor_params_t){ .film_stock = DT_FILMSTOCK_NB,
                                                                 .Dmin = { 1.0f, 1.0f, 1.0f, 0.0f},
                                                                 .wb_high = { 1.0f, 1.0f, 1.0f, 0.0f },
                                                                 .wb_low = { 1.0f, 1.0f, 1.0f, 0.0f },
                                                                 .D_max = 2.2f,
                                                                 .offset = -0.05f,
                                                                 .gamma = 5.0f,
                                                                 .soft_clip = 0.75f,
                                                                 .exposure = 1.f,
                                                                 .black = 0.0755f,
                                                                 .checker = 0};


  dt_gui_presets_add_generic(_("black and white film"), self->op,
                             self->version(), &tmq, sizeof(tmq), 1, DEVELOP_BLEND_CS_RGB_DISPLAY);
}

void init_global(dt_iop_module_so_t *module)
{
  dt_iop_negadoctor_global_data_t *gd
      = (dt_iop_negadoctor_global_data_t *)malloc(sizeof(dt_iop_negadoctor_global_data_t));

  module->data = gd;
  const int program = 30; // negadoctor.cl, from programs.conf
  gd->kernel_negadoctor = dt_opencl_create_kernel(program, "negadoctor");
}

void cleanup_global(dt_iop_module_so_t *module)
{
  dt_iop_negadoctor_global_data_t *gd = module->data;
  dt_opencl_free_kernel(gd->kernel_negadoctor);
  free(module->data);
  module->data = NULL;
}

void init_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  piece->data = g_malloc0(sizeof(dt_iop_negadoctor_data_t));
  piece->data_size = sizeof(dt_iop_negadoctor_data_t);
}

void cleanup_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  g_free(piece->data);
  piece->data = NULL;
}


/* Global GUI stuff */

static void setup_color_variables(dt_iop_negadoctor_gui_data_t *const g, const gint state)
{
  gtk_widget_set_visible(g->Dmin_G, state);
  gtk_widget_set_visible(g->Dmin_B, state);
}


static void toggle_stock_controls(dt_iop_module_t *const self)
{
  dt_iop_negadoctor_gui_data_t *const g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  const dt_iop_negadoctor_params_t *const p = (dt_iop_negadoctor_params_t *)self->params;

  if(p->film_stock == DT_FILMSTOCK_NB)
  {
    // Hide color controls
    setup_color_variables(g, FALSE);
    dt_bauhaus_widget_set_label(g->Dmin_R, NULL, N_("D min"));
  }
  else if(p->film_stock == DT_FILMSTOCK_COLOR)
  {
    // Show color controls
    setup_color_variables(g, TRUE);
    dt_bauhaus_widget_set_label(g->Dmin_R, NULL, N_("D min red component"));
  }
  else
  {
    // We shouldn't be there
    fprintf(stderr, "negadoctor film stock: undefined behaviour\n");
  }
}


static void Dmin_picker_update(dt_iop_module_t *self)
{
  dt_iop_negadoctor_gui_data_t *const g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  const dt_iop_negadoctor_params_t *const p = (dt_iop_negadoctor_params_t *)self->params;

  GdkRGBA color;
  color.alpha = 1.0f;

  if(p->film_stock == DT_FILMSTOCK_COLOR)
  {
    color.red = p->Dmin[0];
    color.green = p->Dmin[1];
    color.blue = p->Dmin[2];
  }
  else if(p->film_stock == DT_FILMSTOCK_NB)
  {
    color.red = color.green = color.blue = p->Dmin[0];
  }

  gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(g->Dmin_picker), &color);
}

static void Dmin_picker_callback(GtkColorButton *widget, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  dt_iop_color_picker_reset(self, TRUE);

  GdkRGBA c;
  gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(widget), &c);
  p->Dmin[0] = c.red;
  p->Dmin[1] = c.green;
  p->Dmin[2] = c.blue;

  ++darktable.gui->reset;
  dt_bauhaus_slider_set(g->Dmin_R, p->Dmin[0]);
  dt_bauhaus_slider_set(g->Dmin_G, p->Dmin[1]);
  dt_bauhaus_slider_set(g->Dmin_B, p->Dmin[2]);
  --darktable.gui->reset;

  Dmin_picker_update(self);
  dt_iop_color_picker_reset(self, TRUE);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void WB_low_picker_update(dt_iop_module_t *self)
{
  dt_iop_negadoctor_gui_data_t *const g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  const dt_iop_negadoctor_params_t *const p = (dt_iop_negadoctor_params_t *)self->params;

  GdkRGBA color;
  color.alpha = 1.0f;

  dt_aligned_pixel_t WB_low_invert;
  for(size_t c = 0; c < 3; ++c) WB_low_invert[c] = 2.0f - p->wb_low[c];
  const float WB_low_max = v_maxf(WB_low_invert);
  for(size_t c = 0; c < 3; ++c) WB_low_invert[c] /= WB_low_max;

  color.red = WB_low_invert[0];
  color.green = WB_low_invert[1];
  color.blue = WB_low_invert[2];

  gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(g->WB_low_picker), &color);
}

static void WB_low_picker_callback(GtkColorButton *widget, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  dt_iop_color_picker_reset(self, TRUE);

  GdkRGBA c;
  gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(widget), &c);

  dt_aligned_pixel_t RGB = { 2.0f - c.red, 2.0f - c.green, 2.0f - c.blue };

  float RGB_min = v_minf(RGB);
  for(size_t k = 0; k < 3; k++) p->wb_low[k] = RGB[k] / RGB_min;

  ++darktable.gui->reset;
  dt_bauhaus_slider_set(g->wb_low_R, p->wb_low[0]);
  dt_bauhaus_slider_set(g->wb_low_G, p->wb_low[1]);
  dt_bauhaus_slider_set(g->wb_low_B, p->wb_low[2]);
  --darktable.gui->reset;

  WB_low_picker_update(self);
  dt_iop_color_picker_reset(self, TRUE);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}


static void WB_high_picker_update(dt_iop_module_t *self)
{
  dt_iop_negadoctor_gui_data_t *const g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  const dt_iop_negadoctor_params_t *const p = (dt_iop_negadoctor_params_t *)self->params;

  GdkRGBA color;
  color.alpha = 1.0f;

  dt_aligned_pixel_t WB_high_invert;
  for(size_t c = 0; c < 3; ++c) WB_high_invert[c] = 2.0f - p->wb_high[c];
  const float WB_high_max = v_maxf(WB_high_invert);
  for(size_t c = 0; c < 3; ++c) WB_high_invert[c] /= WB_high_max;

  color.red = WB_high_invert[0];
  color.green = WB_high_invert[1];
  color.blue = WB_high_invert[2];

  gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(g->WB_high_picker), &color);
}

static void WB_high_picker_callback(GtkColorButton *widget, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  dt_iop_color_picker_reset(self, TRUE);

  GdkRGBA c;
  gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(widget), &c);

  dt_aligned_pixel_t RGB = { 2.0f - c.red, 2.0f - c.green, 2.0f - c.blue };
  float RGB_min = v_minf(RGB);
  for(size_t k = 0; k < 3; k++) p->wb_high[k] = RGB[k] / RGB_min;

  ++darktable.gui->reset;
  dt_bauhaus_slider_set(g->wb_high_R, p->wb_high[0]);
  dt_bauhaus_slider_set(g->wb_high_G, p->wb_high[1]);
  dt_bauhaus_slider_set(g->wb_high_B, p->wb_high[2]);
  --darktable.gui->reset;

  WB_high_picker_update(self);
  dt_iop_color_picker_reset(self, TRUE);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void Wb_low_norm_callback(GtkWidget *widget, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;

  dt_iop_negadoctor_gui_data_t *const g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;


  const float WB_low_max = v_maxf(p->wb_low);
  for(size_t c = 0; c < 3; ++c)
    p->wb_low[c] -= WB_low_max;


  ++darktable.gui->reset;
  dt_bauhaus_slider_set(g->wb_low_R, p->wb_low[0]);
  dt_bauhaus_slider_set(g->wb_low_G, p->wb_low[1]);
  dt_bauhaus_slider_set(g->wb_low_B, p->wb_low[2]);
  --darktable.gui->reset;

  WB_low_picker_update(self);
  dt_control_queue_redraw_widget(self->widget);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void Wb_high_norm_callback(GtkWidget *widget, dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;

  dt_iop_negadoctor_gui_data_t *const g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  const float WB_high_min = v_minf(p->wb_high);
  for(size_t c = 0; c < 3; ++c)
    p->wb_high[c] /= WB_high_min;


  ++darktable.gui->reset;
  dt_bauhaus_slider_set(g->wb_high_R, p->wb_high[0]);
  dt_bauhaus_slider_set(g->wb_high_G, p->wb_high[1]);
  dt_bauhaus_slider_set(g->wb_high_B, p->wb_high[2]);
  --darktable.gui->reset;

  WB_low_picker_update(self);
  dt_control_queue_redraw_widget(self->widget);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

/* Color pickers auto-tuners */

// measure Dmin from the film edges first
static void apply_auto_Dmin(dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  for(int k = 0; k < 4; k++) p->Dmin[k] = self->picked_color[k];

  ++darktable.gui->reset;
  dt_bauhaus_slider_set(g->Dmin_R, p->Dmin[0]);
  dt_bauhaus_slider_set(g->Dmin_G, p->Dmin[1]);
  dt_bauhaus_slider_set(g->Dmin_B, p->Dmin[2]);
  --darktable.gui->reset;

  Dmin_picker_update(self);
  dt_control_queue_redraw_widget(self->widget);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

// from Dmin, find out the range of density values of the film and compute Dmax
static void apply_auto_Dmax(dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  dt_aligned_pixel_t RGB;
  for(int c = 0; c < 3; c++)
  {
    RGB[c] = log10f(p->Dmin[c] / fmaxf(self->picked_color_min[c], THRESHOLD));
  }

  // Take the max(RGB) for safety. Big values unclip whites
  p->D_max = v_maxf(RGB);

  ++darktable.gui->reset;
  dt_bauhaus_slider_set(g->D_max, p->D_max);
  --darktable.gui->reset;

  dt_control_queue_redraw_widget(self->widget);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

// from Dmax, compute the offset so the range of density is rescaled between [0; 1]
static void apply_auto_offset(dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  dt_aligned_pixel_t RGB;
  for(int c = 0; c < 3; c++)
    RGB[c] = log10f(p->Dmin[c] / fmaxf(self->picked_color_max[c], THRESHOLD)) / p->D_max;

  // Take the min(RGB) for safety. Negative values unclip blacks
  p->offset = v_minf(RGB);

  ++darktable.gui->reset;
  dt_bauhaus_slider_set(g->offset, p->offset);
  --darktable.gui->reset;

  dt_control_queue_redraw_widget(self->widget);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

// from Dmax and offset, compute the white balance correction as multipliers of the offset
// such that offset x wb[c] make black monochrome
static void apply_auto_WB_low(dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  dt_aligned_pixel_t RGB_min;
  for(int c = 0; c < 3; c++)
    RGB_min[c] = p->offset - log10f(p->Dmin[c] / fmaxf(self->picked_color[c], THRESHOLD)) / p->D_max;

  const float RGB_v_min = v_minf(RGB_min); // warning: can be negative
  for(int c = 0; c < 3; c++) p->wb_low[c] = (RGB_min[c] - RGB_v_min) / -0.3f;

  fprintf(stdout,"NEG offset: %f\n", p->offset); 
  fprintf(stdout,"NEG LOW: %f %f %f\n",RGB_min[0], RGB_min[1], RGB_min[2]);
  fprintf(stdout,"NEG LOW res: %f %f %f\n",RGB_min[0]- RGB_v_min, RGB_min[1]- RGB_v_min, RGB_min[2]- RGB_v_min);

  ++darktable.gui->reset;
  dt_bauhaus_slider_set(g->wb_low_R, p->wb_low[0]);
  dt_bauhaus_slider_set(g->wb_low_G, p->wb_low[1]);
  dt_bauhaus_slider_set(g->wb_low_B, p->wb_low[2]);
  --darktable.gui->reset;

  WB_low_picker_update(self);
  dt_control_queue_redraw_widget(self->widget);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

// from Dmax, offset and white balance multipliers, compute the white balance of the illuminant as multipliers of 1/Dmax
// such that WB[c] / Dmax make white monochrome
static void apply_auto_WB_high(dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  dt_aligned_pixel_t RGB_min;
  for(int c = 0; c < 3; c++)
    RGB_min[c] = fabsf(-1.0f / ((p->offset + p->wb_low[c] * -0.3f) - log10f(p->Dmin[c] / fmaxf(self->picked_color[c], THRESHOLD)) / p->D_max));

  const float RGB_v_min = v_minf(RGB_min); // warning : must be positive
  for(int c = 0; c < 3; c++) p->wb_high[c] = RGB_min[c] / RGB_v_min;

  ++darktable.gui->reset;
  dt_bauhaus_slider_set(g->wb_high_R, p->wb_high[0]);
  dt_bauhaus_slider_set(g->wb_high_G, p->wb_high[1]);
  dt_bauhaus_slider_set(g->wb_high_B, p->wb_high[2]);
  --darktable.gui->reset;

  WB_high_picker_update(self);
  dt_control_queue_redraw_widget(self->widget);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

// from Dmax, offset and both white balances, compute the print black adjustment
// such that the printed values range from 0 to + infinity
static void apply_auto_black(dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  dt_aligned_pixel_t RGB;
  for(int c = 0; c < 3; c++)
  {
    RGB[c] = -log10f(p->Dmin[c] / fmaxf(self->picked_color_max[c], THRESHOLD));
    RGB[c] *= p->wb_high[c] / p->D_max;
    RGB[c] += p->wb_low[c] * (p->offset + p->wb_low[c] * -0.3f);
    RGB[c] = 0.1f - (1.0f - fast_exp10f(RGB[c])); // actually, remap between -3.32 EV and infinity for safety because gamma comes later
  }
  p->black = v_maxf(RGB);

  ++darktable.gui->reset;
  dt_bauhaus_slider_set(g->black, p->black);
  --darktable.gui->reset;

  dt_control_queue_redraw_widget(self->widget);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

// from Dmax, offset, both white balances, and printblack, compute the print exposure adjustment as a scaling factor
// such that the printed values range from 0 to 1
static void apply_auto_exposure(dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  dt_aligned_pixel_t RGB;
  for(int c = 0; c < 3; c++)
  {
    RGB[c] = -log10f(p->Dmin[c] / fmaxf(self->picked_color_min[c], THRESHOLD));
    RGB[c] *= p->wb_high[c] / p->D_max;
    RGB[c] += (p->offset + p->wb_low[c] * -0.3f);
    RGB[c] = 0.96f / (1.0f - fast_exp10f(RGB[c]) + p->black); // actually, remap in [0; 0.96] for safety
  }
  p->exposure = v_minf(RGB);

  ++darktable.gui->reset;
  dt_bauhaus_slider_set(g->exposure, log2f(p->exposure));
  --darktable.gui->reset;

  dt_control_queue_redraw_widget(self->widget);
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void _do_rgb_densities(dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  const dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  const dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  float offset = p->D_film_max[0] / p->D[0];
  float thing[3] = { 0.f };
  float average[2][3] = {{ 0.f }};

  for(int k = 0; k < 9; k++) 
  {
    thing[k % 3] = p->D[k] * offset;
    gtk_label_set_text(GTK_LABEL(g->density_info[k]), g_strdup_printf("%.2f", thing[k % 3]));

    if(k % 3 == 2)
    {
      average[0][k / 3] = v_minf(thing);
      average[1][k / 3] = v_maxf(thing);
    }
  }

  for(int a = 0; a < 3; a++)
  {
    float averaged = (average[0][a] + average[1][a]) / 2;
    gtk_label_set_text(GTK_LABEL(g->density_info[a + 9]), g_strdup_printf("%.2f", averaged));
  }
}

static void apply_dmin_range_picker(dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  ++darktable.gui->reset;
  for(int c = 0; c < 3; c++) 
  {
    p->D[c] = log10f(1.f / self->picked_color[c]);
    fprintf(stdout,"NEG Dmin[%i]: %.2f\n", c, p->D[c]);
  }
  --darktable.gui->reset;

  _do_rgb_densities(self);
}

static void apply_dmax_range_picker(dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  ++darktable.gui->reset;
  fprintf(stdout,"NEG sizeof D: %i\n", sizeof(p->D)/sizeof(float));

  for(int c = 0; c < 3; c++) 
  {
    p->D[c + 3] = log10f(1.f / self->picked_color[c]);
    fprintf(stdout,"NEG D[%i]: %.2f\n",c + 3, p->D[c + 3]);
  }
  --darktable.gui->reset;

  _do_rgb_densities(self);
}

static void density_picker(dt_iop_module_t *self)
{
  if(darktable.gui->reset) return;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  ++darktable.gui->reset;
  for(int c = 0; c < 3; c++)
    p->D[c + 6] = log10f(1.f / self->picked_color[c]);
  --darktable.gui->reset;

  _do_rgb_densities(self);
  
}

void color_picker_apply(dt_iop_module_t *self, GtkWidget *picker, dt_dev_pixelpipe_iop_t *piece)
{
  if(darktable.gui->reset || picker == NULL) return;
  const dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;

  if     (picker == g->Dmin_sampler)
    apply_auto_Dmin(self);
  else if(picker == g->WB_high_sampler)
    apply_auto_WB_high(self);
  else if(picker == g->offset)
    apply_auto_offset(self);
  else if(picker == g->D_max)
    apply_auto_Dmax(self);
  else if(picker == g->WB_low_sampler)
    apply_auto_WB_low(self);
  else if(picker == g->exposure)
    apply_auto_exposure(self);
  else if(picker == g->black)
    apply_auto_black(self);
  else if(picker == g->D_sampler)
    density_picker(self);
  else if(picker == g->Dmin_range_sampler)
    apply_dmin_range_picker(self);
  else if(picker == g->Dmax_range_sampler)
    apply_dmax_range_picker(self);
  else
    fprintf(stderr, "[negadoctor] unknown color picker\n");
}

//copied from import.c
static GtkWidget * _attach_aligned_grid_item(GtkWidget *grid, const int row, const int column,
                                      const char *label, const GtkAlign align, const gboolean fixed_width,
                                      const gboolean full_width)
{
  GtkWidget *w = gtk_label_new(label);
  if(fixed_width)
    gtk_label_set_max_width_chars(GTK_LABEL(w), 25);

  gtk_label_set_ellipsize(GTK_LABEL(w), PANGO_ELLIPSIZE_END);
  gtk_grid_attach(GTK_GRID(grid), w, column, row, full_width ? 2 : 1, 1);
  gtk_label_set_xalign(GTK_LABEL(w), align);
  gtk_widget_set_halign(w, align);
  gtk_label_set_line_wrap(GTK_LABEL(w), TRUE);
  return w;
}

static inline void update_bounding_box(dt_iop_negadoctor_gui_data_t *g,
                                       const float x_increment, const float y_increment)
{
  // update box nodes
  for(size_t k = 0; k < 4; k++)
  {
    if(g->active_node[k])
    {
      g->box[k].x += x_increment;
      g->box[k].y += y_increment;
    }
  }

  // update the homography
  get_homography(g->ideal_box, g->box, g->homography);
  get_homography(g->box, g->ideal_box, g->inverse_homography);
}

static inline void init_bounding_box(dt_iop_negadoctor_gui_data_t *g, const float width)
{
  if(!g->checker_ready)
  {
    // top left
    g->box[0].x = g->box[0].y = 10.;

    // top right
    g->box[1].x = width - 10.f;
    g->box[1].y = g->box[0].y;

    // bottom right
    g->box[2].x = g->box[1].x;
    g->box[2].y = (width - 10.f) * g->checker->ratio;

    // bottom left
    g->box[3].x = g->box[0].x;
    g->box[3].y = g->box[2].y;

    g->checker_ready = TRUE;
  }

  g->center_box.x = 0.5f;
  g->center_box.y = 0.5f;

  g->ideal_box[0].x = 0.f;
  g->ideal_box[0].y = 0.f;
  g->ideal_box[1].x = 1.f;
  g->ideal_box[1].y = 0.f;
  g->ideal_box[2].x = 1.f;
  g->ideal_box[2].y = 1.f;
  g->ideal_box[3].x = 0.f;
  g->ideal_box[3].y = 1.f;

  update_bounding_box(g, 0.f, 0.f);
}

static void _start_profiling(dt_iop_module_t *self, const gint page_num)
{
  if(darktable.gui->reset) return;
  
  fprintf(stdout,"START PROF: BEGIN\n");
  if(!self || self == NULL)
  {
    fprintf(stdout,"START PROF: self is NULL\n");
    return;
  }
  fprintf(stdout,"START PROF: self exists\n");
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;
  dt_iop_request_focus(self);
  fprintf(stdout,"START PROF: page %i\n", page_num);
  fprintf(stdout,"START PROF: checker: %i\n", p->checker);
  
  if(p->checker && page_num == 3)
  {
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(self->off), TRUE);

    const dt_develop_t *dev = self->dev;
    const float wd = (float) dev->preview_pipe->backbuf_width;
    const float ht = (float) dev->preview_pipe->backbuf_height;
    if(wd == 0.f || ht == 0.f) return;

    g->is_profiling_started = TRUE;
    fprintf(stdout,"START PROF: profiling will start\n");
    // init bounding box
    dt_iop_gui_enter_critical_section(self);
    init_bounding_box(g, wd);
    dt_iop_gui_leave_critical_section(self);
  }
  else
  {
    fprintf(stdout,"START PROF: Deactivate profiling\n");
    g->is_profiling_started = FALSE;
  }

  fprintf(stdout,"START PROF: redraw center\n");
  dt_control_queue_redraw_center();
  fprintf(stdout,"START PROF: redraw center OK\n");
  
}


static void start_profiling_callback(GtkNotebook* notebook, GtkWidget* page, guint page_num, gpointer user_data)
{
  if(darktable.gui->reset) return;
  dt_iop_module_t *self = (dt_iop_module_t *) user_data;
  if(!self || self == NULL)
  {
    fprintf(stdout,"callback: self is NULL\n");
    return;
  }
  fprintf(stdout,"callback: notebook page\n");
  _start_profiling(self, page_num);
  fprintf(stdout,"callback: notebook page end\n");
}

int mouse_moved(struct dt_iop_module_t *self, double x, double y, double pressure, int which)
{
  if(!self->enabled) return 0;

  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  if(g == NULL || !g->is_profiling_started) return 0;
  if(g->box[0].x == -1.0f || g->box[1].y == -1.0f) return 0;

  dt_develop_t *dev = self->dev;
  const float wd = dev->preview_pipe->backbuf_width;
  const float ht = dev->preview_pipe->backbuf_height;
  if(wd == 0.f || ht == 0.f) return 0;

  float pzx, pzy;
  dt_dev_get_pointer_zoom_pos(dev, x, y, &pzx, &pzy);
  pzx += 0.5f;
  pzy += 0.5f;
  pzx *= wd;
  pzy *= ht;

  // if dragging and dropping, don't update active nodes,
  // just update cursor coordinates then redraw
  // this ensure smooth updates
  if(g->drag_drop)
  {
    dt_iop_gui_enter_critical_section(self);
    g->click_end.x = pzx;
    g->click_end.y = pzy;

    update_bounding_box(g, g->click_end.x - g->click_start.x, g->click_end.y - g->click_start.y);

    g->click_start.x = pzx;
    g->click_start.y = pzy;
    dt_iop_gui_leave_critical_section(self);

    dt_control_queue_redraw_center();
    return 1;
  }

  // Find out if we are close to a node
  dt_iop_gui_enter_critical_section(self);
  g->is_cursor_close = FALSE;

  for(size_t k = 0; k < 4; k++)
  {
    if(hypotf(pzx - g->box[k].x, pzy - g->box[k].y) < 15.f)
    {
      g->active_node[k] = TRUE;
      g->is_cursor_close = TRUE;
    }
    else
      g->active_node[k] = FALSE;
  }
  dt_iop_gui_leave_critical_section(self);

  // if cursor is close from a node, remove the system pointer arrow to prevent hiding the spot behind it
  if(g->is_cursor_close)
  {
    dt_control_change_cursor(GDK_BLANK_CURSOR);
  }
  else
  {
    // fall back to default cursor
    GdkCursor *const cursor = gdk_cursor_new_from_name(gdk_display_get_default(), "default");
    gdk_window_set_cursor(gtk_widget_get_window(dt_ui_main_window(darktable.gui->ui)), cursor);
    g_object_unref(cursor);
  }

  dt_control_queue_redraw_center();

  return 1;
}

int button_pressed(struct dt_iop_module_t *self, double x, double y, double pressure, int which, int type,
                   uint32_t state)
{
  if(!self->enabled) return 0;

  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  if(g == NULL || !g->is_profiling_started) return 0;

  dt_develop_t *dev = self->dev;
  const float wd = dev->preview_pipe->backbuf_width;
  const float ht = dev->preview_pipe->backbuf_height;
  if(wd == 0.f || ht == 0.f) return 0;

  // double click : reset the perspective correction
  if(type == GDK_DOUBLE_BUTTON_PRESS)
  {
    dt_iop_gui_enter_critical_section(self);
    g->checker_ready = FALSE;
    //g->profile_ready = FALSE;
    init_bounding_box(g, wd);
    dt_iop_gui_leave_critical_section(self);

    dt_control_queue_redraw_center();
    return 1;
  }

  // bounded box not inited, abort
  if(g->box[0].x == -1.0f || g->box[1].y == -1.0f) return 0;

  // cursor is not on a node, abort
  if(!g->is_cursor_close) return 0;

  float pzx, pzy;
  dt_dev_get_pointer_zoom_pos(dev, x, y, &pzx, &pzy);
  pzx += 0.5f;
  pzy += 0.5f;
  pzx *= wd;
  pzy *= ht;

  dt_iop_gui_enter_critical_section(self);
  g->drag_drop = TRUE;
  g->click_start.x = pzx;
  g->click_start.y = pzy;
  dt_iop_gui_leave_critical_section(self);

  dt_control_queue_redraw_center();

  return 1;
}

int button_released(struct dt_iop_module_t *self, double x, double y, int which, uint32_t state)
{
  if(!self->enabled) return 0;

  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  if(g == NULL || !g->is_profiling_started) return 0;
  if(g->box[0].x == -1.0f || g->box[1].y == -1.0f) return 0;
  if(!g->is_cursor_close || !g->drag_drop) return 0;

  dt_develop_t *dev = self->dev;
  const float wd = dev->preview_pipe->backbuf_width;
  const float ht = dev->preview_pipe->backbuf_height;
  if(wd == 0.f || ht == 0.f) return 0;

  float pzx, pzy;
  dt_dev_get_pointer_zoom_pos(dev, x, y, &pzx, &pzy);
  pzx += 0.5f;
  pzy += 0.5f;
  pzx *= wd;
  pzy *= ht;

  dt_iop_gui_enter_critical_section(self);
  g->drag_drop = FALSE;
  g->click_end.x = pzx;
  g->click_end.y = pzy;
  update_bounding_box(g, g->click_end.x - g->click_start.x, g->click_end.y - g->click_start.y);
  dt_iop_gui_leave_critical_section(self);

  dt_control_queue_redraw_center();

  return 1;
}

void gui_post_expose(struct dt_iop_module_t *self, cairo_t *cr, int32_t width, int32_t height,
                     int32_t pointerx, int32_t pointery)
{
  const dt_iop_order_iccprofile_info_t *const work_profile = dt_ioppr_get_pipe_output_profile_info(self->dev->pipe);
  if(work_profile == NULL) return;

  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  if(!g->is_profiling_started) return;

  // Rescale and shift Cairo drawing coordinates
  dt_develop_t *dev = self->dev;
  const float wd = dev->preview_pipe->backbuf_width;
  const float ht = dev->preview_pipe->backbuf_height;
  if(wd == 0.f || ht == 0.f) return;

  const float zoom_y = dt_control_get_dev_zoom_y();
  const float zoom_x = dt_control_get_dev_zoom_x();
  const dt_dev_zoom_t zoom = dt_control_get_dev_zoom();
  const int closeup = dt_control_get_dev_closeup();
  const float zoom_scale = dt_dev_get_zoom_scale(dev, zoom, 1<<closeup, 1);
  cairo_translate(cr, width / 2.0, height / 2.0);
  cairo_scale(cr, zoom_scale, zoom_scale);
  cairo_translate(cr, -.5f * wd - zoom_x * wd, -.5f * ht - zoom_y * ht);

  cairo_set_line_width(cr, 2.0 / zoom_scale);
  const double origin = 9. / zoom_scale;
  const double destination = 18. / zoom_scale;

  for(size_t k = 0; k < 4; k++)
  {
    if(g->active_node[k])
    {
      // draw cross hair
      cairo_set_source_rgba(cr, 1., 1., 1., 1.);

      cairo_move_to(cr, g->box[k].x - origin, g->box[k].y);
      cairo_line_to(cr, g->box[k].x - destination, g->box[k].y);

      cairo_move_to(cr, g->box[k].x + origin, g->box[k].y);
      cairo_line_to(cr, g->box[k].x + destination, g->box[k].y);

      cairo_move_to(cr, g->box[k].x, g->box[k].y - origin);
      cairo_line_to(cr, g->box[k].x, g->box[k].y - destination);

      cairo_move_to(cr, g->box[k].x, g->box[k].y + origin);
      cairo_line_to(cr, g->box[k].x, g->box[k].y + destination);

      cairo_stroke(cr);
    }

    // draw outline circle
    cairo_set_source_rgba(cr, 1., 1., 1., 1.);
    cairo_arc(cr, g->box[k].x, g->box[k].y, 8. / zoom_scale, 0, 2. * M_PI);
    cairo_stroke(cr);

    // draw black dot
    cairo_set_source_rgba(cr, 0., 0., 0., 1.);
    cairo_arc(cr, g->box[k].x, g->box[k].y, 1.5 / zoom_scale, 0, 2. * M_PI);
    cairo_fill(cr);
  }

  // draw symmetry axes
  cairo_set_line_width(cr, 1.5 / zoom_scale);
  cairo_set_source_rgba(cr, 1., 1., 1., 1.);
  const point_t top_ideal = { 0.5f, 1.f };
  const point_t top = apply_homography(top_ideal, g->homography);
  const point_t bottom_ideal = { 0.5f, 0.f };
  const point_t bottom = apply_homography(bottom_ideal, g->homography);
  cairo_move_to(cr, top.x, top.y);
  cairo_line_to(cr, bottom.x, bottom.y);
  cairo_stroke(cr);

  const point_t left_ideal = { 0.f, 0.5f };
  const point_t left = apply_homography(left_ideal, g->homography);
  const point_t right_ideal = { 1.f, 0.5f };
  const point_t right = apply_homography(right_ideal, g->homography);
  cairo_move_to(cr, left.x, left.y);
  cairo_line_to(cr, right.x, right.y);
  cairo_stroke(cr);

  /* For debug : display center of the image and center of the ideal target
  point_t new_target_center = apply_homography(target_center, g->homography);
  cairo_set_source_rgba(cr, 1., 1., 1., 1.);
  cairo_arc(cr, new_target_center.x, new_target_center.y, 7., 0, 2. * M_PI);
  cairo_stroke(cr);

  cairo_set_source_rgba(cr, 0., 1., 1., 1.);
  cairo_arc(cr, 0.5 * wd, 0.5 * ht, 7., 0, 2. * M_PI);
  cairo_stroke(cr);
  */

  const float radius_x = g->checker->radius * hypotf(1.f, g->checker->ratio) * g->safety_margin;
  const float radius_y = radius_x / g->checker->ratio;

  for(size_t k = 0; k < g->checker->patches; k++)
  {
    // center of the patch in the ideal reference
    const point_t center = { g->checker->values[k].x, g->checker->values[k].y };

    // corners of the patch in the ideal reference
    const point_t corners[4] = { {center.x - radius_x, center.y - radius_y},
                                 {center.x + radius_x, center.y - radius_y},
                                 {center.x + radius_x, center.y + radius_y},
                                 {center.x - radius_x, center.y + radius_y} };

    // apply patch coordinates transform depending on perspective
    const point_t new_center = apply_homography(center, g->homography);
    // apply_homography_scaling gives a scaling of areas. we need to scale the
    // radius of the center circle so take a square root.
    const float scaling = sqrtf(apply_homography_scaling(center, g->homography));
    point_t new_corners[4];
    for(size_t c = 0; c < 4; c++) new_corners[c] = apply_homography(corners[c], g->homography);

    cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
    cairo_set_source_rgba(cr, 0., 0., 0., 1.);
    cairo_move_to(cr, new_corners[0].x, new_corners[0].y);
    cairo_line_to(cr, new_corners[1].x, new_corners[1].y);
    cairo_line_to(cr, new_corners[2].x, new_corners[2].y);
    cairo_line_to(cr, new_corners[3].x, new_corners[3].y);
    cairo_line_to(cr, new_corners[0].x, new_corners[0].y);

    /*if(g->delta_E_in)
    {
      // draw delta E feedback
      if(g->delta_E_in[k] > 2.3f)
      {
        // one diagonal if delta E > 3
        cairo_move_to(cr, new_corners[0].x, new_corners[0].y);
        cairo_line_to(cr, new_corners[2].x, new_corners[2].y);
      }
      if(g->delta_E_in[k] > 4.6f)
      {
        // the other diagonal if delta E > 6
        cairo_move_to(cr, new_corners[1].x, new_corners[1].y);
        cairo_line_to(cr, new_corners[3].x, new_corners[3].y);
      }
    }
*/
    cairo_set_line_width(cr, 5.0 / zoom_scale);
    cairo_stroke_preserve(cr);
    cairo_set_line_width(cr, 2.0 / zoom_scale);
    cairo_set_source_rgba(cr, 1., 1., 1., 1.);
    cairo_stroke(cr);

    cairo_set_line_cap(cr, CAIRO_LINE_CAP_BUTT);

    dt_aligned_pixel_t RGB;
    dt_ioppr_lab_to_rgb_matrix(g->checker->values[k].Lab, RGB, work_profile->matrix_out_transposed, work_profile->lut_out,
                               work_profile->unbounded_coeffs_out, work_profile->lutsize,
                               work_profile->nonlinearlut);

    cairo_set_source_rgba(cr, RGB[0], RGB[1], RGB[2], 1.);
    cairo_arc(cr, new_center.x, new_center.y, 0.25 * (radius_x + radius_y) * scaling, 0, 2. * M_PI);
    cairo_fill(cr);
  }
}

static void checker_changed_callback(GtkWidget *widget, gpointer user_data)
{
  if(darktable.gui->reset) return;
  dt_iop_module_t *self = (dt_iop_module_t *)user_data;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;

  p->checker = dt_bauhaus_combobox_get(widget); //get choosen step wedge 
  fprintf(stdout, "CHECKER number: %i\n", p->checker);
  if(p->checker == 0)
  {
    fprintf(stdout, "CHECKER: NONE\n");
    gtk_widget_set_visible(g->safety, FALSE);
    gtk_widget_set_visible(g->button_commit, FALSE);
    g->is_profiling_started = FALSE;
    g->checker_ready = FALSE;
    dt_control_queue_redraw_center();
  }
  else
  {
    g->checker = dt_get_negadoctor_step_wedge(p->checker - 1); // use selected step wedge
    gtk_widget_set_visible(g->safety, TRUE);
    gtk_widget_set_visible(g->button_commit, TRUE);
  
    // now we have a step wedge
    fprintf(stdout, "CHECKER: %s\n", g->checker->name);
    _start_profiling(self, gtk_notebook_get_current_page(g->notebook));
  }
  //dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void safety_changed_callback(GtkWidget *widget, gpointer user_data)
{
  if(darktable.gui->reset) return;
  dt_iop_module_t *self = (dt_iop_module_t *)user_data;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;

  dt_iop_gui_enter_critical_section(self);
  g->safety_margin = dt_bauhaus_slider_get(widget);
  dt_iop_gui_leave_critical_section(self);

  dt_conf_set_float("darkroom/modules/negadoctor/safety", g->safety_margin);
  dt_control_queue_redraw_center();
}

static void run_validation_callback(GtkWidget *widget, GdkEventButton *event, gpointer user_data)
{
  if(darktable.gui->reset) return;
  dt_iop_module_t *self = (dt_iop_module_t *)user_data;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;

  dt_iop_gui_enter_critical_section(self);
  g->run_validation = TRUE;
  dt_iop_gui_leave_critical_section(self);

  dt_dev_invalidate_preview(self->dev);
  dt_dev_refresh_ui_images(self->dev);
}


void gui_init(dt_iop_module_t *self)
{
  dt_iop_negadoctor_gui_data_t *g = IOP_GUI_ALLOC(negadoctor);

 // Init the color checker UI
  for(size_t k = 0; k < 4; k++)
  {
    g->box[k].x = g->box[k].y = -1.;
    g->active_node[k] = FALSE;
  }
  g->is_cursor_close = FALSE;
  g->drag_drop = FALSE;
  g->checker_ready = FALSE;

  static dt_action_def_t notebook_def = { };
  g->notebook = dt_ui_notebook_new(&notebook_def);
  dt_action_define_iop(self, NULL, N_("page"), GTK_WIDGET(g->notebook), &notebook_def);
  g_signal_connect(G_OBJECT(g->notebook), "switch-page", G_CALLBACK(start_profiling_callback), self);

  // Page FILM PROPERTIES
  GtkWidget *page1 = self->widget = dt_ui_notebook_page(g->notebook, N_("film properties"), NULL);

  // Dmin

  gtk_box_pack_start(GTK_BOX(page1), dt_ui_section_label_new(_("color of the film base")), FALSE, FALSE, 0);

  GtkWidget *row1 = GTK_WIDGET(gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0));

  g->Dmin_picker = gtk_color_button_new();
  gtk_color_chooser_set_use_alpha(GTK_COLOR_CHOOSER(g->Dmin_picker), FALSE);
  gtk_color_button_set_title(GTK_COLOR_BUTTON(g->Dmin_picker), _("select color of film material from a swatch"));
  gtk_box_pack_start(GTK_BOX(row1), GTK_WIDGET(g->Dmin_picker), TRUE, TRUE, 0);
  g_signal_connect(G_OBJECT(g->Dmin_picker), "color-set", G_CALLBACK(Dmin_picker_callback), self);

  g->Dmin_sampler = dt_color_picker_new(self, DT_COLOR_PICKER_AREA, row1);
  gtk_widget_set_tooltip_text(g->Dmin_sampler , _("pick color of film material from image"));

  gtk_box_pack_start(GTK_BOX(page1), GTK_WIDGET(row1), FALSE, FALSE, 0);

  g->Dmin_R = dt_bauhaus_slider_from_params(self, "Dmin[0]");
  dt_bauhaus_slider_set_digits(g->Dmin_R, 4);
  dt_bauhaus_slider_set_format(g->Dmin_R, "%");
  dt_bauhaus_slider_set_factor(g->Dmin_R, 100);
  dt_bauhaus_widget_set_label(g->Dmin_R, NULL, N_("D min red component"));
  gtk_widget_set_tooltip_text(g->Dmin_R, _("adjust the color and shade of the film transparent base.\n"
                                           "this value depends on the film material, \n"
                                           "the chemical fog produced while developing the film,\n"
                                           "and the scanner white balance."));

  g->Dmin_G = dt_bauhaus_slider_from_params(self, "Dmin[1]");
  dt_bauhaus_slider_set_digits(g->Dmin_G, 4);
  dt_bauhaus_slider_set_format(g->Dmin_G, "%");
  dt_bauhaus_slider_set_factor(g->Dmin_G, 100);
  dt_bauhaus_widget_set_label(g->Dmin_G, NULL, N_("D min green component"));
  gtk_widget_set_tooltip_text(g->Dmin_G, _("adjust the color and shade of the film transparent base.\n"
                                           "this value depends on the film material, \n"
                                           "the chemical fog produced while developing the film,\n"
                                           "and the scanner white balance."));

  g->Dmin_B = dt_bauhaus_slider_from_params(self, "Dmin[2]");
  dt_bauhaus_slider_set_digits(g->Dmin_B, 4);
  dt_bauhaus_slider_set_format(g->Dmin_B, "%");
  dt_bauhaus_slider_set_factor(g->Dmin_B, 100);
  dt_bauhaus_widget_set_label(g->Dmin_B, NULL, N_("D min blue component"));
  gtk_widget_set_tooltip_text(g->Dmin_B, _("adjust the color and shade of the film transparent base.\n"
                                           "this value depends on the film material, \n"
                                           "the chemical fog produced while developing the film,\n"
                                           "and the scanner white balance."));

  // D max and scanner bias

  gtk_box_pack_start(GTK_BOX(page1), dt_ui_section_label_new(_("dynamic range of the film")), FALSE, FALSE, 0);

  g->D_max = dt_color_picker_new(self, DT_COLOR_PICKER_AREA, dt_bauhaus_slider_from_params(self, "D_max"));
  dt_bauhaus_slider_set_format(g->D_max, " dB");
  gtk_widget_set_tooltip_text(g->D_max, _("maximum density of the film, corresponding to white after inversion.\n"
                                          "this value depends on the film specifications, the developing process,\n"
                                          "the dynamic range of the scene and the scanner exposure settings."));

  gtk_box_pack_start(GTK_BOX(page1), dt_ui_section_label_new(_("scanner exposure settings")), FALSE, FALSE, 0);

  g->offset = dt_color_picker_new(self, DT_COLOR_PICKER_AREA, dt_bauhaus_slider_from_params(self, "offset"));
  dt_bauhaus_slider_set_format(g->offset, " dB");
  gtk_widget_set_tooltip_text(g->offset, _("correct the exposure of the scanner, for all RGB channels,\n"
                                           "before the inversion, so blacks are neither clipped or too pale."));

  // Page CORRECTIONS
  GtkWidget *page2 = self->widget = dt_ui_notebook_page(g->notebook, N_("corrections"), NULL);

  // WB shadows
  gtk_box_pack_start(GTK_BOX(page2), dt_ui_section_label_new(_("shadows color cast")), FALSE, FALSE, 0);

  GtkWidget *row3 = GTK_WIDGET(gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0));

  g->WB_low_picker = gtk_color_button_new();
  gtk_color_chooser_set_use_alpha(GTK_COLOR_CHOOSER(g->WB_low_picker), FALSE);
  gtk_color_button_set_title(GTK_COLOR_BUTTON(g->WB_low_picker), _("select color of shadows from a swatch"));
  gtk_box_pack_start(GTK_BOX(row3), GTK_WIDGET(g->WB_low_picker), TRUE, TRUE, 0);
  g_signal_connect(G_OBJECT(g->WB_low_picker), "color-set", G_CALLBACK(WB_low_picker_callback), self);

  g->WB_low_norm = dt_action_button_new((dt_lib_module_t *)self, N_("normalize"), Wb_low_norm_callback, self, _("normalize shadows white balance settings"), 0, 0);
  gtk_box_pack_start(GTK_BOX(row3), GTK_WIDGET(g->WB_low_norm), FALSE, FALSE, 0);

  g->WB_low_sampler = dt_color_picker_new(self, DT_COLOR_PICKER_AREA, row3);
  gtk_widget_set_tooltip_text(g->WB_low_sampler, _("pick shadows color from image"));

  gtk_box_pack_start(GTK_BOX(page2), GTK_WIDGET(row3), FALSE, FALSE, 0);

  g->wb_low_R = dt_bauhaus_slider_from_params(self, "wb_low[0]");
  dt_bauhaus_widget_set_label(g->wb_low_R, NULL, N_("shadows red offset"));
  gtk_widget_set_tooltip_text(g->wb_low_R, _("correct the color cast in shadows so blacks are\n"
                                             "truly achromatic. Setting this value before\n"
                                             "the highlights illuminant white balance will help\n"
                                             "recovering the global white balance in difficult cases."));

  g->wb_low_G = dt_bauhaus_slider_from_params(self, "wb_low[1]");
  dt_bauhaus_widget_set_label(g->wb_low_G, NULL, N_("shadows green offset"));
  gtk_widget_set_tooltip_text(g->wb_low_G, _("correct the color cast in shadows so blacks are\n"
                                             "truly achromatic. Setting this value before\n"
                                             "the highlights illuminant white balance will help\n"
                                             "recovering the global white balance in difficult cases."));

  g->wb_low_B = dt_bauhaus_slider_from_params(self, "wb_low[2]");
  dt_bauhaus_widget_set_label(g->wb_low_B, NULL, N_("shadows blue offset"));
  gtk_widget_set_tooltip_text(g->wb_low_B, _("correct the color cast in shadows so blacks are\n"
                                             "truly achromatic. Setting this value before\n"
                                             "the highlights illuminant white balance will help\n"
                                             "recovering the global white balance in difficult cases."));

  // WB highlights
  gtk_box_pack_start(GTK_BOX(page2), dt_ui_section_label_new(_("highlights white balance")), FALSE, FALSE, 0);

  GtkWidget *row2 = GTK_WIDGET(gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0));

  g->WB_high_picker = gtk_color_button_new();
  gtk_color_chooser_set_use_alpha(GTK_COLOR_CHOOSER(g->WB_high_picker), FALSE);
  gtk_color_button_set_title(GTK_COLOR_BUTTON(g->WB_high_picker), _("select color of illuminant from a swatch"));
  gtk_box_pack_start(GTK_BOX(row2), GTK_WIDGET(g->WB_high_picker), TRUE, TRUE, 0);
  g_signal_connect(G_OBJECT(g->WB_high_picker), "color-set", G_CALLBACK(WB_high_picker_callback), self);

  g->WB_high_norm = dt_action_button_new((dt_lib_module_t *)self, N_("normalize"), Wb_high_norm_callback, self, _("normalize highlight white balance settings"), 0, 0);
  gtk_box_pack_start(GTK_BOX(row2), GTK_WIDGET(g->WB_high_norm), FALSE, FALSE, 0);

  g->WB_high_sampler = dt_color_picker_new(self, DT_COLOR_PICKER_AREA, row2);
  gtk_widget_set_tooltip_text(g->WB_high_sampler , _("pick illuminant color from image"));

  gtk_box_pack_start(GTK_BOX(page2), GTK_WIDGET(row2), FALSE, FALSE, 0);

  g->wb_high_R = dt_bauhaus_slider_from_params(self, "wb_high[0]");
  dt_bauhaus_widget_set_label(g->wb_high_R, NULL, N_("illuminant red gain"));
  gtk_widget_set_tooltip_text(g->wb_high_R, _("correct the color of the illuminant so whites are\n"
                                              "truly achromatic. Setting this value after\n"
                                              "the shadows color cast will help\n"
                                              "recovering the global white balance in difficult cases."));

  g->wb_high_G = dt_bauhaus_slider_from_params(self, "wb_high[1]");
  dt_bauhaus_widget_set_label(g->wb_high_G, NULL, N_("illuminant green gain"));
  gtk_widget_set_tooltip_text(g->wb_high_G, _("correct the color of the illuminant so whites are\n"
                                              "truly achromatic. Setting this value after\n"
                                              "the shadows color cast will help\n"
                                              "recovering the global white balance in difficult cases."));

  g->wb_high_B = dt_bauhaus_slider_from_params(self, "wb_high[2]");
  dt_bauhaus_widget_set_label(g->wb_high_B, NULL, N_("illuminant blue gain"));
  gtk_widget_set_tooltip_text(g->wb_high_B, _("correct the color of the illuminant so whites are\n"
                                              "truly achromatic. Setting this value after\n"
                                              "the shadows color cast will help\n"
                                              "recovering the global white balance in difficult cases."));

  // Page PRINT PROPERTIES
  GtkWidget *page3 = self->widget = dt_ui_notebook_page(g->notebook, N_("print properties"), NULL);

  // print corrections
  gtk_box_pack_start(GTK_BOX(page3), dt_ui_section_label_new(_("virtual paper properties")), FALSE, FALSE, 0);

  g->black = dt_color_picker_new(self, DT_COLOR_PICKER_AREA, dt_bauhaus_slider_from_params(self, "black"));
  dt_bauhaus_slider_set_digits(g->black, 4);
  dt_bauhaus_slider_set_factor(g->black, 100);
  dt_bauhaus_slider_set_format(g->black, "%");
  gtk_widget_set_tooltip_text(g->black, _("correct the density of black after the inversion,\n"
                                          "to adjust the global contrast while avoiding clipping shadows."));

  g->gamma = dt_bauhaus_slider_from_params(self, "gamma");
  dt_bauhaus_widget_set_label(g->gamma, NULL, N_("paper grade (gamma)"));
  gtk_widget_set_tooltip_text(g->gamma, _("select the grade of the virtual paper, which is actually\n"
                                          "equivalent to applying a gamma. it compensates the film D max\n"
                                          "and recovers the contrast. use a high grade for high D max."));

  g->soft_clip = dt_bauhaus_slider_from_params(self, "soft_clip");
  dt_bauhaus_slider_set_factor(g->soft_clip, 100);
  dt_bauhaus_slider_set_digits(g->soft_clip, 4);
  dt_bauhaus_slider_set_format(g->soft_clip, "%");
  gtk_widget_set_tooltip_text(g->soft_clip, _("gradually compress specular highlights past this value\n"
                                              "to avoid clipping while pushing the exposure for mid-tones.\n"
                                              "this somewhat reproduces the behaviour of matte paper."));

  gtk_box_pack_start(GTK_BOX(page3), dt_ui_section_label_new(_("virtual print emulation")), FALSE, FALSE, 0);

  g->exposure = dt_color_picker_new(self, DT_COLOR_PICKER_AREA, dt_bauhaus_slider_from_params(self, "exposure"));
  dt_bauhaus_slider_set_soft_min(g->exposure, -1.0);
  dt_bauhaus_slider_set_soft_max(g->exposure, 1.0);
  dt_bauhaus_slider_set_default(g->exposure, 0.0);
  dt_bauhaus_slider_set_format(g->exposure, _(" EV"));
  gtk_widget_set_tooltip_text(g->exposure, _("correct the printing exposure after inversion to adjust\n"
                                             "the global contrast and avoid clipping highlights."));

  // Densitometer

  GtkWidget *page4 = self->widget = dt_ui_notebook_page(g->notebook, N_("Density"), NULL);

  gtk_box_pack_start(GTK_BOX(page4), dt_ui_section_label_new(_("Step wedge")), FALSE, FALSE, 0);

  GtkWidget *g_stepwedge = self->widget = GTK_WIDGET(gtk_box_new(GTK_ORIENTATION_VERTICAL, 0));

  DT_BAUHAUS_COMBOBOX_NEW_FULL(g->checkers_list, self, NULL, N_("Use a step wedge"),
                                _("choose the vendor and the type of your chart"),
                                0, checker_changed_callback, self,
                                N_("None"),
                                N_("Stouffer 21 Steps Wedge"));
  gtk_box_pack_start(GTK_BOX(g_stepwedge), GTK_WIDGET(g->checkers_list), TRUE, FALSE, 0);

  g->safety = dt_bauhaus_slider_new_with_range_and_feedback(self, 0., 1., 0, 0.5, 3, TRUE);
  dt_bauhaus_widget_set_label(g->safety, NULL, N_("patch scale"));
  gtk_widget_set_tooltip_text(g->safety, _("reduce the radius of the patches to select the more or less central part.\n"
                                           "useful when the perspective correction is sloppy or\n"
                                           "the patches frame cast a shadows on the edges of the patch." ));
  g_signal_connect(G_OBJECT(g->safety), "value-changed", G_CALLBACK(safety_changed_callback), self);
  gtk_box_pack_start(GTK_BOX(g_stepwedge), GTK_WIDGET(g->safety), TRUE, FALSE, 0);

  GtkWidget *toolbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, DT_BAUHAUS_SPACE);

  g->button_commit = dtgtk_button_new(dtgtk_cairo_paint_check_mark, 0, NULL);
  g_signal_connect(G_OBJECT(g->button_commit), "button-press-event", G_CALLBACK(run_validation_callback), (gpointer)self);
  gtk_widget_set_tooltip_text(g->button_commit, _("Read density value"));
  gtk_box_pack_end(GTK_BOX(toolbar), GTK_WIDGET(g->button_commit), FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(g_stepwedge), GTK_WIDGET(toolbar), FALSE, FALSE, 0);

 
  gtk_box_pack_start(GTK_BOX(page4), GTK_WIDGET(g_stepwedge), FALSE, FALSE, 0);

  gtk_box_pack_start(GTK_BOX(page4), dt_ui_section_label_new(_("Densities")), FALSE, FALSE, 0);

  GtkWidget *g_densito = self->widget = GTK_WIDGET(gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0));

  //g->Dmin_range_sampler = dt_color_picker_new(self, DT_COLOR_PICKER_AREA, GTK_WIDGET(range_row));
  g->Dmin_range_sampler = dt_color_picker_new(self, DT_COLOR_PICKER_AREA, NULL);
  gtk_widget_set_tooltip_text(g->Dmin_range_sampler , _("Dmin Range"));

  g->Dmax_range_sampler = dt_color_picker_new(self, DT_COLOR_PICKER_AREA, NULL);
  gtk_widget_set_tooltip_text(g->Dmin_range_sampler , _("Dmax Range"));

  GtkWidget *range = gtk_grid_new();
  gtk_grid_set_column_spacing(GTK_GRID(range), 10);
  gtk_grid_set_row_spacing(GTK_GRID(range), 1);
  _attach_aligned_grid_item(range, 1, 0, _("Red:"), GTK_ALIGN_END, FALSE, FALSE);
  _attach_aligned_grid_item(range, 2, 0, _("Green:"), GTK_ALIGN_END, FALSE, FALSE);
  _attach_aligned_grid_item(range, 3, 0, _("Blue:"), GTK_ALIGN_END, FALSE, FALSE);
  _attach_aligned_grid_item(range, 5, 0, _("Average:"), GTK_ALIGN_END, FALSE, FALSE);

  _attach_aligned_grid_item(range, 0, 1, _("Dmin"), GTK_ALIGN_START, FALSE, FALSE);
  g->density_info[RANGE_MIN_RED] = _attach_aligned_grid_item(range, 1, 1, "", GTK_ALIGN_START, TRUE, FALSE);
  g->density_info[RANGE_MIN_GREEN] = _attach_aligned_grid_item(range, 2, 1, "", GTK_ALIGN_START, TRUE, FALSE);
  g->density_info[RANGE_MIN_BLUE] = _attach_aligned_grid_item(range, 3, 1, "", GTK_ALIGN_START, TRUE, FALSE);
  gtk_grid_attach(GTK_GRID(range),  g->Dmax_range_sampler, 1, 4, 1, 1);
  g->density_info[RANGE_MIN_AVERAGE] = _attach_aligned_grid_item(range, 5, 1, "", GTK_ALIGN_START, TRUE, FALSE);

  _attach_aligned_grid_item(range, 0, 2, _("Dmax"), GTK_ALIGN_START, FALSE, FALSE);
  g->density_info[RANGE_MAX_RED] = _attach_aligned_grid_item(range, 1, 2, "", GTK_ALIGN_START, TRUE, FALSE);
  g->density_info[RANGE_MAX_GREEN] = _attach_aligned_grid_item(range, 2, 2, "", GTK_ALIGN_START, TRUE, FALSE);
  g->density_info[RANGE_MAX_BLUE] = _attach_aligned_grid_item(range, 3, 2, "", GTK_ALIGN_START, TRUE, FALSE);
  gtk_grid_attach(GTK_GRID(range),  g->Dmin_range_sampler, 2, 4, 1, 1);
  g->density_info[RANGE_MAX_AVERAGE] = _attach_aligned_grid_item(range, 5, 2, "", GTK_ALIGN_START, TRUE, FALSE);

  gtk_box_pack_start(GTK_BOX(g_densito), GTK_WIDGET(range), TRUE, FALSE, 0);

  
  g->D_sampler = dt_color_picker_new(self, DT_COLOR_PICKER_POINT_AREA, NULL);
  gtk_widget_set_tooltip_text(g->D_sampler , _("Read the density value"));

  GtkWidget *spot = gtk_grid_new();
  gtk_grid_set_column_spacing(GTK_GRID(spot), 10);
  gtk_grid_set_row_spacing(GTK_GRID(spot), 1);

  _attach_aligned_grid_item(spot, 1, 0, _("Red:"), GTK_ALIGN_END, FALSE, FALSE);
  _attach_aligned_grid_item(spot, 2, 0, _("Green:"), GTK_ALIGN_END, FALSE, FALSE);
  _attach_aligned_grid_item(spot, 3, 0, _("Blue:"), GTK_ALIGN_END, FALSE, FALSE);
  _attach_aligned_grid_item(spot, 5, 0, _("Average:"), GTK_ALIGN_END, FALSE, FALSE);

  _attach_aligned_grid_item(spot, 0, 1, _("Spot"), GTK_ALIGN_START, FALSE, FALSE);
  g->density_info[SPOT_RED] = _attach_aligned_grid_item(spot, 1, 1, "", GTK_ALIGN_START, TRUE, FALSE);
  g->density_info[SPOT_GREEN] = _attach_aligned_grid_item(spot, 2, 1, "", GTK_ALIGN_START, TRUE, FALSE);
  g->density_info[SPOT_BLUE] = _attach_aligned_grid_item(spot, 3, 1, "", GTK_ALIGN_START, TRUE, FALSE);
  gtk_grid_attach(GTK_GRID(spot), g->D_sampler, 1, 4, 1, 1);
  g->density_info[SPOT_AVERAGE] = _attach_aligned_grid_item(spot, 5, 1, "", GTK_ALIGN_START, TRUE, FALSE);

  gtk_box_pack_start(GTK_BOX(g_densito), GTK_WIDGET(spot), TRUE, FALSE, 0);

  gtk_box_pack_start(GTK_BOX(page4), GTK_WIDGET(g_densito), FALSE, FALSE, 0);

  self->widget = GTK_WIDGET(page4);
  g->D_film_max_R = dt_bauhaus_slider_from_params(self, "D_film_max[0]");
  dt_bauhaus_slider_set_digits(g->D_film_max_R, 2);
  dt_bauhaus_slider_set_format(g->D_film_max_R, " D");
  dt_bauhaus_widget_set_label(g->D_film_max_R, NULL, N_("Film Dmax red"));
  gtk_widget_set_tooltip_text(g->D_film_max_R, _("Set the Dmax of the film to calibrate values.\n"));

  //for(int k = 0; k < 12; k++) gtk_label_set_text(GTK_LABEL(g->density_info[k]), "n/a");
    
  // start building top level widget
  self->widget = gtk_box_new(GTK_ORIENTATION_VERTICAL, DT_BAUHAUS_SPACE);

  // Film emulsion
  g->film_stock = dt_bauhaus_combobox_from_params(self, "film_stock");
  gtk_widget_set_tooltip_text(g->film_stock, _("toggle on or off the color controls"));

  gtk_box_pack_start(GTK_BOX(self->widget), GTK_WIDGET(g->notebook), FALSE, FALSE, 0);
}


void gui_changed(dt_iop_module_t *self, GtkWidget *w, void *previous)
{
  dt_iop_negadoctor_params_t *p = (dt_iop_negadoctor_params_t *)self->params;
  dt_iop_negadoctor_gui_data_t *g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  if(!w || w == g->film_stock)
  {
    toggle_stock_controls(self);
    Dmin_picker_update(self);
  }
  else if(w == g->Dmin_R && p->film_stock == DT_FILMSTOCK_NB)
  {
    dt_bauhaus_slider_set(g->Dmin_G, p->Dmin[0]);
    dt_bauhaus_slider_set(g->Dmin_B, p->Dmin[0]);
  }
  else if(w == g->Dmin_R || w == g->Dmin_G || w == g->Dmin_B)
  {
    Dmin_picker_update(self);
  }
  else if(w == g->exposure)
  {
    p->exposure = powf(2.0f, p->exposure);
  }

  if(!w || w == g->wb_high_R || w == g->wb_high_G || w == g->wb_high_B)
  {
    WB_high_picker_update(self);
  }

  if(!w || w == g->wb_low_R || w == g->wb_low_G || w == g->wb_low_B)
  {
    WB_low_picker_update(self);
  }

  if(!w || w == g->D_film_max_R)
  {
    _do_rgb_densities(self);
  }    
}


void gui_update(dt_iop_module_t *const self)
{
  // let gui slider match current parameters:
  dt_iop_negadoctor_gui_data_t *const g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  const dt_iop_negadoctor_params_t *const p = (dt_iop_negadoctor_params_t *)self->params;

  dt_iop_color_picker_reset(self, TRUE);
  dt_bauhaus_slider_set(g->exposure, log2f(p->exposure));     // warning: GUI is in EV

  // Update custom stuff
  dt_iop_gui_enter_critical_section(self);
  dt_bauhaus_combobox_set(g->checkers_list, p->checker);
  if(p->checker)
  {
    gtk_widget_set_visible(g->safety, TRUE);
    gtk_widget_set_visible(g->button_commit, TRUE);
  }
  else
  {
    gtk_widget_set_visible(g->safety, FALSE);
    gtk_widget_set_visible(g->button_commit, FALSE);
  }
  g->checker = dt_get_negadoctor_step_wedge(p->checker);
  g->safety_margin = dt_conf_get_float("darkroom/modules/negadoctor/safety");
  dt_bauhaus_slider_set(g->safety, g->safety_margin);
  dt_iop_gui_leave_critical_section(self);

  gui_changed(self, NULL, NULL);
}

void gui_reset(dt_iop_module_t *self)
{
  dt_iop_negadoctor_gui_data_t *const g = (dt_iop_negadoctor_gui_data_t *)self->gui_data;
  dt_iop_color_picker_reset(self, TRUE);
  g->is_profiling_started = FALSE;
  gui_changed(self, NULL, NULL);
}
// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
