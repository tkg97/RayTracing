// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <GLFW/glfw3.h>
GLFWwindow* window;

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
using namespace glm;

#include "common/shader.hpp"
#include "common/texture.hpp"
#include "common/objloader.hpp"

int render(const std::vector<float> &rayData)
{
	// Initialise GLFW
	if( !glfwInit() )
	{
		fprintf( stderr, "Failed to initialize GLFW\n" );
		getchar();
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow( 1024, 768, "Visualizer", NULL, NULL);
	if( window == NULL ){
		fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
		getchar();
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		glfwTerminate();
		return -1;
	}

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS); 

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	GLuint programID = LoadShaders( "shaders/VertexShading.shader", "shaders/FragmentShading.shader" );
	GLuint programIDline = LoadShaders("shaders/LineVertexShader.shader", "shaders/LineFragmentShader.shader");

	// Get a handle for our "MVP" uniform
	GLuint MatrixID = glGetUniformLocation(programID, "MVP");
	GLuint MatrixIDline = glGetUniformLocation(programID, "MVP");
	GLuint ViewMatrixID = glGetUniformLocation(programID, "V");
	GLuint ModelMatrixID = glGetUniformLocation(programID, "M");
	// Get a handle for our "LightPosition" uniform
	GLuint LightID = glGetUniformLocation(programID, "LightPosition_worldspace");

	// Load the texture
	GLuint TextureSphere = loadBMP_custom("inputFiles/Opengl/texture_red.bmp");
	GLuint TexturePlane = loadBMP_custom("inputFiles/Opengl/texture_grey.bmp");
	
	// Get a handle for our "myTextureSampler" uniform
	GLuint TextureID  = glGetUniformLocation(programID, "myTextureSampler");

	// Read our .obj file
	std::vector<glm::vec3> verticesSphere, verticesPlane;
	std::vector<glm::vec2> uvsSphere, uvsPlane;
	std::vector<glm::vec3> normalsSphere, normalsPlane;
	bool res = loadOBJ("inputFiles/Opengl/sphere.obj", verticesSphere, uvsSphere, normalsSphere);

	if (!res) {
		std::cout << "Sphere object file couldn't be loaded" << std::endl;
		exit(0);
	}

	res = loadOBJ("inputFiles/Opengl/plane.obj", verticesPlane, uvsPlane, normalsPlane);

	if (!res) {
		std::cout << "Plane object file couldn't be loaded" << std::endl;
	}

	// Load it into a VBO Sphere

	GLuint vertexbufferSphere;
	glGenBuffers(1, &vertexbufferSphere);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbufferSphere);
	glBufferData(GL_ARRAY_BUFFER, verticesSphere.size() * sizeof(glm::vec3), &verticesSphere[0], GL_STATIC_DRAW);

	GLuint uvbufferSphere;
	glGenBuffers(1, &uvbufferSphere);
	glBindBuffer(GL_ARRAY_BUFFER, uvbufferSphere);
	glBufferData(GL_ARRAY_BUFFER, uvsSphere.size() * sizeof(glm::vec2), &uvsSphere[0], GL_STATIC_DRAW);

	GLuint normalbufferSphere;
	glGenBuffers(1, &normalbufferSphere);
	glBindBuffer(GL_ARRAY_BUFFER, normalbufferSphere);
	glBufferData(GL_ARRAY_BUFFER, normalsSphere.size() * sizeof(glm::vec3), &normalsSphere[0], GL_STATIC_DRAW);

	// Load it into a VBO Plane

	GLuint vertexbufferPlane;
	glGenBuffers(1, &vertexbufferPlane);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbufferPlane);
	glBufferData(GL_ARRAY_BUFFER, verticesSphere.size() * sizeof(glm::vec3), &verticesPlane[0], GL_STATIC_DRAW);

	GLuint uvbufferPlane;
	glGenBuffers(1, &uvbufferPlane);
	glBindBuffer(GL_ARRAY_BUFFER, uvbufferPlane);
	glBufferData(GL_ARRAY_BUFFER, uvsPlane.size() * sizeof(glm::vec2), &uvsPlane[0], GL_STATIC_DRAW);

	GLuint normalbufferPlane;
	glGenBuffers(1, &normalbufferPlane);
	glBindBuffer(GL_ARRAY_BUFFER, normalbufferPlane);
	glBufferData(GL_ARRAY_BUFFER, normalsPlane.size() * sizeof(glm::vec3), &normalsPlane[0], GL_STATIC_DRAW);

	// VBO for line
	GLuint vertexbufferLine;
	glGenBuffers(1, &vertexbufferLine);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbufferLine);
	glBufferData(GL_ARRAY_BUFFER, rayData.size() * sizeof(float), &rayData[0], GL_STATIC_DRAW);

	do{

		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Use our shader
		glUseProgram(programID);

		// render sphere first

		glm::mat4 ProjectionMatrix = glm::perspective<float>(glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
		glm::mat4 ViewMatrix = glm::lookAt(glm::vec3(0,0,-20), glm::vec3(0,0,1), glm::vec3(0,1,0));
		
		glm::mat4 ModelMatrixSphere = glm::mat4(1.0);
		ModelMatrixSphere = glm::translate(ModelMatrixSphere, glm::vec3(3.0f, 0.0f, 5.0f));
		glm::mat4 MVPsphere = ProjectionMatrix * ViewMatrix * ModelMatrixSphere;

		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVPsphere[0][0]);
		glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrixSphere[0][0]);
		glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &ViewMatrix[0][0]);

		glm::vec3 lightPos = glm::vec3(4,4,-4);
		glUniform3f(LightID, lightPos.x, lightPos.y, lightPos.z);

		// Bind our texture in Texture Unit 0
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, TextureSphere);
		// Set our "myTextureSampler" sampler to use Texture Unit 0
		glUniform1i(TextureID, 0);

		// 1rst attribute buffer : vertices
		glBindBuffer(GL_ARRAY_BUFFER, vertexbufferSphere);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(
			0,                  // attribute
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// 2nd attribute buffer : UVs
		glBindBuffer(GL_ARRAY_BUFFER, uvbufferSphere);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(
			1,                                // attribute
			2,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		// 3rd attribute buffer : normals
		glBindBuffer(GL_ARRAY_BUFFER, normalbufferSphere);
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(
			2,                                // attribute
			3,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		// Draw the triangles !
		glDrawArrays(GL_TRIANGLES, 0, verticesSphere.size() );


		// now rendering of plane

		glm::mat4 ModelMatrixPlane = glm::mat4(1.0);
		ModelMatrixPlane = glm::translate(ModelMatrixPlane, glm::vec3(-7.0f, 0.0f, 0.0f));
		ModelMatrixPlane = glm::rotate(ModelMatrixPlane, glm::radians(-45.0f), { 0,1,0 });
		ModelMatrixPlane = glm::rotate(ModelMatrixPlane, glm::radians(-90.0f), { 1,0,0 });
		ModelMatrixPlane = glm::scale(ModelMatrixPlane, vec3(0.5f, 0.5f, 0.5f));
		glm::mat4 MVPplane = ProjectionMatrix * ViewMatrix * ModelMatrixPlane;

		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVPplane[0][0]);
		glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrixPlane[0][0]);
		glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &ViewMatrix[0][0]);

		// Bind our texture in Texture Unit 0
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, TexturePlane);
		// Set our "myTextureSampler" sampler to use Texture Unit 0
		glUniform1i(TextureID, 0);

		// 1rst attribute buffer : vertices
		glBindBuffer(GL_ARRAY_BUFFER, vertexbufferPlane);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(
			0,                  // attribute
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// 2nd attribute buffer : UVs
		glBindBuffer(GL_ARRAY_BUFFER, uvbufferPlane);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(
			1,                                // attribute
			2,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		// 3rd attribute buffer : normals
		glBindBuffer(GL_ARRAY_BUFFER, normalbufferPlane);
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(
			2,                                // attribute
			3,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		// Draw the triangles !
		glDrawArrays(GL_TRIANGLES, 0, verticesPlane.size() );

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
		glDisableVertexAttribArray(2);

		// Now render the lines

		glUseProgram(programIDline);

		glm::mat4 ModelMatrixLine = glm::mat4(1.0);
		glm::mat4 MVPline = ProjectionMatrix * ViewMatrix * ModelMatrixLine;

		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixIDline, 1, GL_FALSE, &MVPline[0][0]);

		glBindBuffer(GL_ARRAY_BUFFER, vertexbufferLine);
		glEnableVertexAttribArray(0);

		glVertexAttribPointer(
			0,                                // attribute
			3,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		glDrawArrays(GL_LINES, 0, (rayData.size() / 3));

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();

	} // Check if the ESC key was pressed or the window was closed
	while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
		   glfwWindowShouldClose(window) == 0 );

	// Cleanup VBO and shader
	glDeleteBuffers(1, &vertexbufferSphere);
	glDeleteBuffers(1, &uvbufferSphere);
	glDeleteBuffers(1, &normalbufferSphere);
	glDeleteBuffers(1, &vertexbufferPlane);
	glDeleteBuffers(1, &uvbufferPlane);
	glDeleteBuffers(1, &normalbufferPlane);
	glDeleteBuffers(1, &vertexbufferLine);
	glDeleteProgram(programIDline);
	glDeleteProgram(programID);
	glDeleteTextures(1, &TextureSphere);
	glDeleteTextures(1, &TexturePlane);
	glDeleteVertexArrays(1, &VertexArrayID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}