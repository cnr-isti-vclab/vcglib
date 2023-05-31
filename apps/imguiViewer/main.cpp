

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include <stdio.h>
#include "FileBrowser/ImGuiFileBrowser.h"
#include "triangle_mesh_type.h"
#include <wrap/gui/trackball.h>
#include <wrap/gl/trimesh.h>

// Interface variables
GLFWwindow* window;
vcg::GlTrimesh<MyTriMesh> glWrap;
vcg::Trackball trackball;
imgui_addons::ImGuiFileBrowser file_dialog; 
std::ostringstream oss; // String stream used for the log

// Interface options
ImVec4 clear_color = ImVec4(0.33f, 0.33f, 0.33f, 1.00f);
ImVec4 mesh_color = ImVec4(0.79f, 0.9f, 0.87f, 1.00f);
float line_thick=1.f;
int selected_rend_mode = 5;
bool showMesh=true;
bool showTrackBall=true;
vcg::GLW::DrawMode DMode=vcg::GLW::DMSmooth;
std::string PathMesh="../../meshes/fertility_tri.off";

// Mesh data
MyTriMesh mesh;

void LoadMesh(std::string path)
{
    oss<<"Loading "<<path.c_str()<<std::endl;
    int ret = vcg::tri::io::Importer<MyTriMesh>::Open(mesh,path.c_str());
    vcg::tri::UpdateBounding<MyTriMesh>::Box(mesh);
    vcg::tri::UpdateNormal<MyTriMesh>::PerVertexNormalizedPerFaceNormalized(mesh);

    
    if (ret==0)
        oss<<"Loaded "<<mesh.fn<<" faces "<<mesh.vn<<" vertices"<<std::endl;
    else
        oss<<"Mesh not loaded"<<std::endl;
}

static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    if(ImGui::GetIO().WantCaptureMouse) return;
    
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    float xscale, yscale;
    glfwGetWindowContentScale(window, &xscale, &yscale);
    xpos*=xscale;
    ypos*=yscale;
    ypos=(display_h)-ypos;
    trackball.MouseMove((int)xpos,(int)ypos);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if(ImGui::GetIO().WantCaptureMouse) return;

    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    float xscale, yscale;
    glfwGetWindowContentScale(window, &xscale, &yscale);
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    xpos*=xscale;
    ypos*=yscale;
    ypos=(display_h)-ypos;

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
        trackball.MouseDown((int)xpos,(int)ypos,vcg::Trackball::BUTTON_LEFT);

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
        trackball.MouseUp((int)xpos,(int)ypos,vcg::Trackball::BUTTON_LEFT);

}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    trackball.MouseWheel(yoffset/4);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_PRESS)
        trackball.ButtonDown(vcg::Trackball::KEY_SHIFT);
    if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_RELEASE)
        trackball.ButtonUp(vcg::Trackball::KEY_SHIFT);

    if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_PRESS)
        trackball.ButtonDown(vcg::Trackball::KEY_CTRL);
    if (key == GLFW_KEY_LEFT_CONTROL && action == GLFW_RELEASE)
        trackball.ButtonUp(vcg::Trackball::KEY_CTRL);

    if (key == GLFW_KEY_LEFT_ALT && action == GLFW_PRESS)
        trackball.ButtonDown(vcg::Trackball::KEY_ALT);
    if (key == GLFW_KEY_LEFT_ALT && action == GLFW_RELEASE)
        trackball.ButtonUp(vcg::Trackball::KEY_ALT);
}

// Now inside any function
void showMainMenu()
{
    bool open = false, save = false;
    if(ImGui::BeginMainMenuBar())
    {
        if (ImGui::BeginMenu("File"))
        {
            if (ImGui::MenuItem("Open", NULL))
                open = true;
            if (ImGui::MenuItem("Quit", NULL))
                exit(0);
            
            ImGui::EndMenu();
        }
        ImGui::EndMainMenuBar();
    }
    
    //Remember the name to ImGui::OpenPopup() and showFileDialog() must be same...
    if(open)
        ImGui::OpenPopup("Open File");
    
    if(file_dialog.showFileDialog("Open File", imgui_addons::ImGuiFileBrowser::DialogMode::OPEN, ImVec2(700, 310), ".obj,.ply,.off"))
    {
        LoadMesh(file_dialog.selected_path);
        PathMesh= file_dialog.selected_path;
    }
}

void InitGLFW_Window()
{
    if (!glfwInit())
    {
        std::cout<<"Error Initializing GLFW"<<std::endl;
        exit(0);
    }

    window = glfwCreateWindow(1280, 720, "Dear ImGui GLFW+OpenGL2 example", NULL, NULL);
    if (window == NULL)
    {
        std::cout<<"Error Initializing GLFW Window"<<std::endl;
        exit(0);
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    // Setup window
    glfwSetErrorCallback(glfw_error_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetKeyCallback(window, key_callback);
}

void InitIMGui()
{
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls

    // Setup Dear ImGui style
    // ImGui::StyleColorsClassic();
    ImGui::StyleColorsDark();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();
}

void SetRenderBar()
{
    ImGui::Begin("Simple VCG Mesh Visualizer",NULL,ImGuiWindowFlags_AlwaysAutoResize);       // Create a window called "Hello, world!" and append into it.

    ImGui::Checkbox("Show Mesh", &showMesh);
    ImGui::SameLine();
    ImGui::Checkbox("Show TBall", &showTrackBall);

    ImGui::ColorEdit3("mesh color", (float*)&mesh_color); // Edit 3 floats representing a color
    ImGui::ColorEdit3("clear color", (float*)&clear_color); // Edit 3 floats representing a color
    const char* items[] = { "Box", "Points", "DMWire",
                            "Hidden", "Flat", "Smooth",
                            "FlatWire","Radar"
                          };

    if (ImGui::BeginCombo("Rendering", items[selected_rend_mode]))
    {
        for (int i = 0; i < IM_ARRAYSIZE(items); i++)
        {
            bool is_selected = (selected_rend_mode == i);
            if (ImGui::Selectable(items[i], is_selected))
            {
                selected_rend_mode = i;

                switch(selected_rend_mode)
                {
                case 0  : DMode=vcg::GLW::DMBox;      break;
                case 1  : DMode=vcg::GLW::DMPoints;   break;
                case 2  : DMode=vcg::GLW::DMWire;     break;
                case 3  : DMode=vcg::GLW::DMHidden;   break;
                case 4  : DMode=vcg::GLW::DMFlat;     break;
                case 5  : DMode=vcg::GLW::DMSmooth;   break;
                case 6  : DMode=vcg::GLW::DMFlatWire; break;
                default : DMode=vcg::GLW::DMFlatWire;
                }
            }
            if (is_selected)
            {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }
    ImGui::SliderFloat("thickness", &line_thick, 0.5f, 20.0f);

    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
    ImGui::End();
}

void SetLogWindow()
{
    ImGui::Begin("Log",NULL,ImGuiWindowFlags_AlwaysAutoResize);
    std::string LogPrint= oss.str();
    ImGui::Text("%s",LogPrint.c_str());
    ImGui::End();
}

void SetUserBar()
{
    ImGui::Begin("User",NULL,ImGuiWindowFlags_AlwaysAutoResize);
    if(ImGui::Button("Reload Mesh"))
    {   
        LoadMesh(PathMesh);
        oss<<"Reloaded mesh";        
    }
    
    if(ImGui::Button("Smooth Mesh"))
    {
        vcg::tri::Smooth<MyTriMesh>::VertexCoordLaplacian(mesh, 1, false);
        vcg::tri::UpdateNormal<MyTriMesh>::PerVertexNormalizedPerFaceNormalized(mesh);
        oss<<"Smoothed mesh";        
    }

    ImGui::End();
}


void GLDrawMesh()
{
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(40, (GLdouble)display_w/(GLdouble)display_h, 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,4.f,   0,0,0,   0,1,0);

    trackball.GetView();
    trackball.Apply();
    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    
    // save modelview matrix
    glPushMatrix();
    vcg::glScale(2.0f/mesh.bbox.Diag());
    vcg::glTranslate(-mesh.bbox.Center());
    if (showMesh)
    {
        glColor4f(mesh_color.x * mesh_color.w, mesh_color.y * mesh_color.w, mesh_color.z * mesh_color.w, mesh_color.w);
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glLineWidth(line_thick);
        glWrap.SetHintParamf(vcg::GLW::HNPPointSize,line_thick);
        glWrap.Draw(DMode,vcg::GLW::CMNone,vcg::GLW::TMNone);
        glPopAttrib();
    }
    glPopMatrix();
    if (showTrackBall) trackball.DrawPostApply();
}



int main(int, char**)
{   
    LoadMesh(PathMesh);

    glWrap.m=&mesh;
    
    InitGLFW_Window();

    InitIMGui();
    glewInit();

    // Main loop
    while (!glfwWindowShouldClose(window))
    {
        // Poll and handle events (inputs, window resize, etc.)
        // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
        // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
        // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
        // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
        glfwPollEvents();

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        {
            SetRenderBar();
            SetLogWindow();
            SetUserBar();
            showMainMenu();
        }
        // Rendering
        ImGui::Render();

        GLDrawMesh();

        // If you are using this code with non-legacy OpenGL header/contexts (which you should not, prefer using imgui_impl_opengl3.cpp!!),
        // you may need to backup/reset/restore other state, e.g. for current shader using the commented lines below.
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

        glfwMakeContextCurrent(window);
        glfwSwapBuffers(window);
    }

    // Cleanup
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}