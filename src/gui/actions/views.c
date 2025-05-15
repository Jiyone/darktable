
#include "common/darktable.h"
#include "common/debug.h"
#include "control/conf.h"
#include "control/control.h"
#include "develop/develop.h"
#include "gui/gtk.h"

#include "gui/actions/menu.h"


gboolean views_active_callback(GtkWidget *menu_item)
{
  // The insensitive view is the one whose name matches the menu item label
  const dt_view_t *current_view = dt_view_manager_get_current_view(darktable.view_manager);
  const char *current_label = get_custom_data(menu_item);
  return !g_strcmp0(current_label, current_view->module_name);
}

gboolean views_sensitive_callback(GtkWidget *menu_item)
{
  // The insensitive view is the one whose name matches the menu item label
  const dt_view_t *current_view = dt_view_manager_get_current_view(darktable.view_manager);
  const char *current_label = get_custom_data(menu_item);
  return g_strcmp0(current_label, current_view->module_name);
}

#define MACRO_VIEW(view)                                                                                          \
  static gboolean view_switch_to_##view(GtkAccelGroup *group, GObject *acceleratable, guint keyval,               \
                                        GdkModifierType mods, GtkWidget *widget)                                  \
  {                                                                                                               \
    if(views_sensitive_callback(widget))                                                                          \
      dt_ctl_switch_mode_to(#view);                                                                               \
    return TRUE;                                                                                                  \
}

MACRO_VIEW(lighttable);
MACRO_VIEW(darkroom);
MACRO_VIEW(print);
MACRO_VIEW(slideshow);
MACRO_VIEW(map);

void append_views(GtkWidget **menus, GList **lists, const dt_menus_t index)
{
  for(GList *view_iter = darktable.view_manager->views; view_iter; view_iter = g_list_next(view_iter))
  {
    dt_view_t *view = (dt_view_t *)view_iter->data;
    if(view->flags() & VIEW_FLAGS_HIDDEN) continue;

    void *callback = NULL;
    if(!g_strcmp0(view->module_name, "lighttable"))
      callback = view_switch_to_lighttable;
    else if(!g_strcmp0(view->module_name, "darkroom"))
      callback = view_switch_to_darkroom;
    else if(!g_strcmp0(view->module_name, "print"))
      callback = view_switch_to_print;
    else if(!g_strcmp0(view->module_name, "slideshow"))
      callback = view_switch_to_slideshow;
    else if(!g_strcmp0(view->module_name, "map"))
      callback = view_switch_to_map;

    add_sub_menu_entry(menus, lists, view->name(view), index, view->module_name, callback,
                       NULL, views_active_callback, views_sensitive_callback,
                       !g_strcmp0(view->module_name, "lighttable") ? GDK_KEY_Escape : 0, 0);

    // Darkroom is not handled in global menu since it needs to be opened with an image ID,
    // so we only handle it from filmstrip and lighttable thumbnails.
    // Map and Print are too niche to bother.
  }
}

/* TODO ?
* The current logic is to execute state callbacks (active, sensisitive, check) on each menu activation,
* in the menu.h:update_menu_entries() function.
* This is inexpensive as long as there are not too many items.
* The other approach is to connect menu.h:update_entry() to signals, e.g.

    DT_DEBUG_CONTROL_SIGNAL_CONNECT(darktable.signals, DT_SIGNAL_VIEWMANAGER_VIEW_CHANGED,
                              G_CALLBACK(update_entry), self);
    DT_DEBUG_CONTROL_SIGNAL_CONNECT(darktable.signals, DT_SIGNAL_VIEWMANAGER_VIEW_CANNOT_CHANGE,
                                    G_CALLBACK(_lib_viewswitcher_view_cannot_change_callback), self);

    DT_DEBUG_CONTROL_SIGNAL_DISCONNECT(darktable.signals, G_CALLBACK(update_entry), self);
    DT_DEBUG_CONTROL_SIGNAL_DISCONNECT(darktable.signals, G_CALLBACK(_lib_viewswitcher_view_cannot_change_callback), self);
*
* So the update happens as soon as the signal is emited, only for the relevant menuitems.
*
* To re-evaluate in the future...
*/
