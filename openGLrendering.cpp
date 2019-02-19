// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

#include "openGLrendering.h"
#include "common/shader.hpp"
#include "common/texture.hpp"
#include "common/objloader.hpp"

int render(const std::vector<float>& originalRayData, const std::vector<float>& shadowRayData,
	const std::vector<float>& reflectedRayData, const std::vector<float>& refractedRayData)
{
	//Initialize GLFW;
	if (initializeGLFW() < 0) exit(-1);

	// Initialize GLEW
	if (initializeGLEW() < 0) exit(-1);

	glfwPollEvents();

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	GLuint programID = LoadShaders("shaders/VertexShading.shader", "shaders/FragmentShading.shader");
	GLuint programIDline = LoadShaders("shaders/LineVertexShader.shader", "shaders/LineFragmentShader.shader");

	// Get a handle for our "MVP" uniform
	GLuint MatrixID = glGetUniformLocation(programID, "MVP");
	GLuint MatrixIDline = glGetUniformLocation(programIDline, "MVP");
	GLuint ViewMatrixID = glGetUniformLocation(programID, "V");
	GLuint ModelMatrixID = glGetUniformLocation(programID, "M");

	// Get a handle for our "LightPosition" uniform
	GLuint LightID = glGetUniformLocation(programID, "LightPosition_worldspace");
	GLuint colorID = glGetUniformLocation(programIDline, "in_color");

	// Load the texture
	GLuint TextureSphere = loadBMP_custom("inputFiles/Opengl/texture_red.bmp");
	GLuint TexturePlane = loadBMP_custom("inputFiles/Opengl/texture_grey.bmp");

	// Get a handle for our "myTextureSampler" uniform
	GLuint TextureID = glGetUniformLocation(programID, "myTextureSampler");

	// Read our .obj file
	std::vector<glm::vec3> verticesSphere, verticesPlane;
	std::vector<glm::vec2> uvsSphere, uvsPlane;
	std::vector<glm::vec3> normalsSphere, normalsPlane;
	bool res = loadOBJ("inputFiles/Opengl/sphere.obj", verticesSphere, uvsSphere, normalsSphere);

	if (!res) exit(-1);

	res = loadOBJ("inputFiles/Opengl/plane.obj", verticesPlane, uvsPlane, normalsPlane);

	if (!res) exit(-1);

	// Load it into a VBO Sphere

	GLuint vertexbufferSphere;
	setupBuffer(vertexbufferSphere, (verticesSphere.size() * sizeof(glm::vec3)), (&verticesSphere[0]));

	GLuint uvbufferSphere;
	setupBuffer(uvbufferSphere, (uvsSphere.size() * sizeof(glm::vec2)), (&uvsSphere[0]));

	GLuint normalbufferSphere;
	setupBuffer(normalbufferSphere, (normalsSphere.size() * sizeof(glm::vec3)), (&normalsSphere[0]));

	// Load it into a VBO Plane

	GLuint vertexbufferPlane;
	setupBuffer(vertexbufferPlane, (verticesPlane.size() * sizeof(glm::vec3)), (&verticesPlane[0]));

	GLuint uvbufferPlane;
	setupBuffer(uvbufferPlane, (uvsPlane.size() * sizeof(glm::vec2)), (&uvsPlane[0]));

	GLuint normalbufferPlane;
	setupBuffer(normalbufferPlane, (normalsPlane.size() * sizeof(glm::vec3)), (&normalsPlane[0]));

	// VBO for lines

	GLuint vertexbufferLineOriginal;
	setupBuffer(vertexbufferLineOriginal, (originalRayData.size() * sizeof(float)), (&originalRayData[0]));

	GLuint vertexbufferLineShadow;
	setupBuffer(vertexbufferLineShadow, (shadowRayData.size() * sizeof(float)), (&shadowRayData[0]));

	GLuint vertexbufferLineReflection;
	setupBuffer(vertexbufferLineReflection, (reflectedRayData.size() * sizeof(float)), (&reflectedRayData[0]));

	GLuint vertexbufferLineRefraction;
	setupBuffer(vertexbufferLineRefraction, (refractedRayData.size() * sizeof(float)), (&refractedRayData[0]));

	glm::mat4 ProjectionMatrix = glm::perspective<float>(glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
	glm::mat4 ViewMatrix = glm::lookAt(glm::vec3(0, 20, 0), glm::vec3(0, -1, 0), glm::vec3(-1, 0, 0));

	glm::mat4 ModelMatrixLine = glm::mat4(1.0);
	glm::mat4 MVPline = ProjectionMatrix * ViewMatrix * ModelMatrixLine;

	glm::mat4 ModelMatrixSphere = glm::mat4(1.0);
	ModelMatrixSphere = glm::scale(ModelMatrixSphere, glm::vec3(0.5f, 0.5f, 0.5f));
	glm::mat4 MVPsphere = ProjectionMatrix * ViewMatrix * ModelMatrixSphere;

	glm::mat4 ModelMatrixPlane = glm::mat4(1.0);
	ModelMatrixPlane = glm::translate(ModelMatrixPlane, glm::vec3(-10.0f, 0.0f, 0.0f));
	ModelMatrixPlane = glm::rotate(ModelMatrixPlane, glm::radians(-45.0f), { 0,1,0 });
	ModelMatrixPlane = glm::rotate(ModelMatrixPlane, glm::radians(-90.0f), { 1,0,0 });
	ModelMatrixPlane = glm::scale(ModelMatrixPlane, vec3(0.5f, 0.5f, 0.5f));
	glm::mat4 MVPplane = ProjectionMatrix * ViewMatrix * ModelMatrixPlane;

	glm::vec3 lightPos = glm::vec3(0, 0, -8);

	do {

		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Now render the lines

		glUseProgram(programIDline);


		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixIDline, 1, GL_FALSE, &MVPline[0][0]);

		glUniform3f(colorID, 1.0f, 0.0f, 0.0f);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbufferLineOriginal);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		glDrawArrays(GL_LINES, 0, (originalRayData.size() / 3));
		glDisableVertexAttribArray(0);

		glUniform3f(colorID, 1.0f, 1.0f, 1.0f);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbufferLineShadow);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		glDrawArrays(GL_LINES, 0, (shadowRayData.size() / 3));
		glDisableVertexAttribArray(0);

		glUniform3f(colorID, 0.0f, 1.0f, 0.0f);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbufferLineReflection);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		glDrawArrays(GL_LINES, 0, (reflectedRayData.size() / 3));
		glDisableVertexAttribArray(0);

		glUniform3f(colorID, 0.0f, 0.0f, 1.0f);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbufferLineRefraction);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		glDrawArrays(GL_LINES, 0, (refractedRayData.size() / 3));
		glDisableVertexAttribArray(0);

		glUseProgram(programID);

		// render sphere first


		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVPsphere[0][0]);
		glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrixSphere[0][0]);
		glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &ViewMatrix[0][0]);

		glUniform3f(LightID, lightPos.x, lightPos.y, lightPos.z);

		// Bind our texture in Texture Unit 0
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, TextureSphere);
		// Set our "myTextureSampler" sampler to use Texture Unit 0
		glUniform1i(TextureID, 0);

		// 1rst attribute buffer : vertices
		glBindBuffer(GL_ARRAY_BUFFER, vertexbufferSphere);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

		// 2nd attribute buffer : UVs
		glBindBuffer(GL_ARRAY_BUFFER, uvbufferSphere);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

		// 3rd attribute buffer : normals
		glBindBuffer(GL_ARRAY_BUFFER, normalbufferSphere);
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

		// Draw the triangles !
		glDrawArrays(GL_TRIANGLES, 0, verticesSphere.size());

		// now rendering of plane

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
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

		// 2nd attribute buffer : UVs
		glBindBuffer(GL_ARRAY_BUFFER, uvbufferPlane);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

		// 3rd attribute buffer : normals
		glBindBuffer(GL_ARRAY_BUFFER, normalbufferPlane);
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

		// Draw the triangles !
		glDrawArrays(GL_TRIANGLES, 0, verticesPlane.size());

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
		glDisableVertexAttribArray(2);

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	// Cleanup VBO and shader
	glDeleteBuffers(1, &vertexbufferSphere);
	glDeleteBuffers(1, &uvbufferSphere);
	glDeleteBuffers(1, &normalbufferSphere);
	glDeleteBuffers(1, &vertexbufferPlane);
	glDeleteBuffers(1, &uvbufferPlane);
	glDeleteBuffers(1, &normalbufferPlane);
	glDeleteBuffers(1, &vertexbufferLineOriginal);
	glDeleteBuffers(1, &vertexbufferLineShadow);
	glDeleteBuffers(1, &vertexbufferLineReflection);
	glDeleteBuffers(1, &vertexbufferLineRefraction);
	glDeleteProgram(programIDline);
	glDeleteProgram(programID);
	glDeleteTextures(1, &TextureSphere);
	glDeleteTextures(1, &TexturePlane);
	glDeleteVertexArrays(1, &VertexArrayID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}