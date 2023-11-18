/*
    This file is part of darktable,
    Copyright (C) 2010-2021 darktable developers.

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
// our includes go first:
#include "bauhaus/bauhaus.h"
#include "develop/imageop.h"
#include "develop/imageop_gui.h"
#include "gui/color_picker_proxy.h"
#include "gui/gtk.h"
#include "iop/iop_api.h"
#include "iop/normalize.h"

#include <gtk/gtk.h>
#include <stdlib.h>

// This is an example implementation of an image operation module that does nothing useful.
// It demonstrates how the different functions work together. To build your own module,
// take all of the functions that are mandatory, stripping them of comments.
// Then add only the optional functions that are required to implement the functionality
// you need. Don't copy default implementations (hint: if you don't need to change or add
// anything, you probably don't need the copy). Make sure you choose descriptive names
// for your fields and variables. The ones given here are just examples; rename them.
//
// To have your module compile and appear in darkroom, add it to CMakeLists.txt, with
//  add_iop(normalize "normalize.c")
// and to iop_order.c, in the initialisation of legacy_order & v30_order with:
//  { {XX.0f }, "normalize", 0},

// This is the version of the module's parameters,
// and includes version information about compile-time dt.
// The first released version should be 1.
DT_MODULE_INTROSPECTION(1, dt_iop_normalize_data_t)

// TODO: some build system to support dt-less compilation and translation!

// Enums used in params_t can have $DESCRIPTIONs that will be used to
// automatically populate a combobox with dt_bauhaus_combobox_from_params.
// They are also used in the history changes tooltip.
// Combobox options will be presented in the same order as defined here.
// These numbers must not be changed when a new version is introduced.

typedef struct dt_iop_normalize_data_t
{
  float factor;         // $MIN: 0.01 $MAX: 1.0 $DEFAULT: 1.0 $DESCRIPTION: "factor"
  gboolean correction;   // $DESCRIPTION: "correction"
} dt_iop_normalize_data_t;

typedef struct dt_iop_normalize_gui_data_t
{
  // Whatever you need to make your gui happy and provide access to widgets between gui_init, gui_update etc.
  // Stored in self->gui_data while in darkroom.
  // To permanently store per-user gui configuration settings, you could use dt_conf_set/_get.
  GtkWidget *factor, *correction; // this is needed by gui_update
} dt_iop_normalize_gui_data_t;

typedef struct dt_iop_normalize_global_data_t
{

} dt_iop_normalize_global_data_t;


const char *name()
{
  return _("normalize");
}

int flags()
{
  return IOP_FLAGS_INCLUDE_IN_STYLES | IOP_FLAGS_SUPPORTS_BLENDING | IOP_FLAGS_ALLOW_TILING;
}

int default_group()
{
  return IOP_GROUP_REPAIR;
}

const char **description(struct dt_iop_module_t *self)
{
  return dt_iop_set_description(self, _("correct uneven surfaces, dust, etc... from a base picture"),
                                      _("corrective"),
                                      _("linear, RGB, scene-referred"),
                                      _("linear, RGB or XYZ"),
                                      _("linear, RGB, scene-referred"));
}

int default_colorspace(dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  return IOP_CS_RGB;
}


void commit_params(dt_iop_module_t *self, dt_iop_params_t *p1, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  memcpy(piece->data, p1, self->params_size);
}


void process(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid, void *const ovoid,
             const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
  dt_iop_normalize_data_t *d = (dt_iop_normalize_data_t *)piece->data;
  dt_iop_normalize_gui_data_t *g = (dt_iop_normalize_gui_data_t *)self->gui_data;


  const float *const restrict in = (float *)ivoid;
  float *const restrict out = (float *)ovoid;

  const size_t ch = piece->colors;

  dt_aligned_pixel_t maxi = {0.f};

  // get_maxi(in, roi_out, ch, maxi);

  #ifdef _OPENMP
  #pragma omp parallel for simd reduction(max : maxi) default(none) \
  dt_omp_firstprivate(in, roi_out, ch) \
  schedule(static) collapse(2)
  #endif

  for(size_t k = 0; k < (size_t)roi_out->height * roi_out->width * ch; k += ch)
    for(size_t c = 0; c < 3; c++)
    {
      maxi[c] = fmaxf(in[k + c], maxi[c]);
    }

  printf("[NORMALIZE] maxi : %f %f %f\n", maxi[0], maxi[1], maxi[2]);



#ifdef _OPENMP
#pragma omp parallel for simd default(none) \
  dt_omp_firstprivate(in, out, roi_in, roi_out, ch, maxi, d) \
  schedule(static) collapse(2) //shared(d)
#endif

  for(size_t k = 0; k < (size_t)roi_out->height * roi_out->width * ch; k += ch)
    for(size_t c = 0; c < 3; c++)
    {
      float corrective_pixel;
      corrective_pixel = 1 - (in[k + c] / maxi[c]);

      out[k + c] = (in[k + c] / maxi[c]) / (corrective_pixel * d->factor);
      out[k + c] *=  maxi[c];

      out[k + 3] = in[k + 3];
    }

printf("[NORMALIZE] out : %f %f %f\n", out[0], out[1], out[2]);
}


/*
void init(dt_iop_module_t *module) // Optional init and cleanup
{
  dt_iop_default_init(module);
  // module->hide_enable_button = 1;  // this disable "on/off" button.
}


void init_global(dt_iop_module_so_t *module)
{
  module->data = malloc(sizeof(dt_iop_normalize_global_data_t));
}

void cleanup(dt_iop_module_t *module)
{
  // Releases any memory allocated in init(module)
  // Implement this function explicitly if the module allocates additional memory besides (default_)params.
  // this is rare.
  free(module->params);
  module->params = NULL;
  free(module->default_params);
  module->default_params = NULL;
}
*/

void cleanup_global(dt_iop_module_so_t *module)
{
  free(module->data);
  module->data = NULL;
}

/** Put your local callbacks here, be sure to make them static so they won't be visible outside this file! */
/*
static void extra_callback(GtkWidget *w, dt_iop_module_t *self)
{
  // this is important to avoid cycles!
  if(darktable.gui->reset) return;

  dt_iop_normalize_data_t *p = (dt_iop_normalize_data_t *)self->params;
  dt_iop_normalize_gui_data_t *g = (dt_iop_normalize_gui_data_t *)self->gui_data;

  float extra = dt_bauhaus_slider_get(w);

  // Setting a widget value will trigger a callback that will update params.
  // If this is not desirable (because it might result in a cycle) then use
  // ++darktable.gui->reset;
  dt_bauhaus_slider_set(g->factor, p->factor + extra);
  // and reverse with --darktable.gui->reset;

  // If any params updated directly, not via a callback, then
  // let core know of the changes
  dt_dev_add_history_item(darktable.develop, self, TRUE);
}*/

/** optional gui callbacks. */
void gui_changed(dt_iop_module_t *self, GtkWidget *w, void *previous)
{
  // If defined, this gets called when any of the introspection based widgets
  // (created with dt_bauhaus_..._from_params) are changed.
  // The updated value from the widget is already set in params.
  // any additional side-effects can be achieved here.
  dt_iop_normalize_data_t *p = (dt_iop_normalize_data_t *)self->params;
  dt_iop_normalize_gui_data_t *g = (dt_iop_normalize_gui_data_t *)self->gui_data;

  // Test which widget was changed.
  // If allowing w == NULL, this can be called from gui_update, so that
  // gui configuration adjustments only need to be dealt with once, here.
  /*
  if(!w || w == g->method)
  {
    gtk_widget_set_visible(g->check, p->method == DT_normalize_SECOND);
  }
  */
  // Widget configurations that don't depend any any current params values should
  // go in reload_defaults (if they depend on the image) or gui_init.
}



/** gui setup and update, these are needed. */
void gui_update(dt_iop_module_t *self)
{
  dt_iop_normalize_gui_data_t *g = (dt_iop_normalize_gui_data_t *)self->gui_data;
  dt_iop_normalize_data_t *p = (dt_iop_normalize_data_t *)self->params;

  dt_bauhaus_slider_set(g->factor, p->factor);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(g->correction), p->correction);


  // Any configuration changes to the gui that depend on field values should be done here,
  // or can be done in gui_changed which can then be called from here with widget == NULL.
  gui_changed(self, NULL, NULL);
}



void gui_init(dt_iop_module_t *self)
{
  dt_iop_normalize_gui_data_t *g = IOP_GUI_ALLOC(normalize);

  // If the first widget is created using a _from_params call, self->widget does not have to
  // be explicitly initialised, as a new vertical box will be created automatically.
  self->widget = gtk_box_new(GTK_ORIENTATION_VERTICAL, DT_BAUHAUS_SPACE);

  // If the field name should be used as label too, it does not need a $DESCRIPTION;
  // mark it for translation here using N_()
  //
  // A colorpicker can be attached to a slider, as here, or put standalone in a box.
  // When a color is picked, color_picker_apply is called with either the slider or the
  // button that triggered it.
  g->factor = dt_color_picker_new(self, DT_COLOR_PICKER_AREA,
              dt_bauhaus_slider_from_params(self, N_("factor")));
  // The initial slider range can be reduced from the introspection $MIN - $MAX
  //dt_bauhaus_slider_set_soft_range(g->factor, -1.f, 1.f);
  // The default step is range/100, but can be changed here
  dt_bauhaus_slider_set_step(g->factor, .1);
  dt_bauhaus_slider_set_digits(g->factor, 2);
  // Additional parameters determine how the value will be shown.
  // dt_bauhaus_slider_set_format(g->factor, "%");
  // For a percentage, use factor 100.
  //dt_bauhaus_slider_set_factor(g->factor, -1.0f);
  //dt_bauhaus_slider_set_offset(g->factor, -1.0f);
  // Tooltips explain the otherwise compact interface
  gtk_widget_set_tooltip_text(g->factor, _("adjust factor"));

  // A combobox linked to struct field will be filled with the values and $DESCRIPTIONs
  // in the struct definition, in the same order. The automatic callback will put the
  // enum value, not the position within the combobox list, in the field.

  g->correction = dt_bauhaus_toggle_from_params(self, "correction");

  gtk_box_pack_start(GTK_BOX(self->widget), GTK_WIDGET(g->correction), TRUE, TRUE, 0);
  //g_signal_connect(G_OBJECT(g->correction, "value-changed", G_CALLBACK(extra_callback), self);
}

void gui_cleanup(dt_iop_module_t *self)
{
  // This only needs to be provided if gui_init allocates any memory or resources besides
  // self->widget and gui_data_t. The default function (if an explicit one isn't provided here)
  // takes care of gui_data_t (and gtk destroys the widget anyway). If you override the default,
  // you have to do whatever you have to do, and also call IOP_GUI_FREE to clean up gui_data_t.

  IOP_GUI_FREE;
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on

