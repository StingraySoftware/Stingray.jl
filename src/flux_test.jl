# 1. Import Flux
using Flux

# 2. Define a simple model (a single linear layer)
#    It takes 1 input feature and produces 1 output feature.
model = Dense(1 => 1)

# You can inspect the model's initial parameters (Weight W and bias b)
println("Initial parameters: ", Flux.params(model))

# 3. Generate some synthetic data
#    Let's try to learn the function y = 2x - 1
xtrain = hcat(0.0f0, 1.0f0, 2.0f0, 3.0f0, 4.0f0) # Input features (as a 1xN matrix)
ytrain = hcat(2*xtrain[1,:] .- 1 ...)          # Target outputs (as a 1xN matrix)
data = [(xtrain, ytrain)]                      # Flux typically expects data as a collection of (input, output) tuples

println("\nTraining data (x): ", xtrain)
println("Target data (y): ", ytrain)

# 4. Define a loss function
#    Mean Squared Error is common for regression
loss(x, y) = Flux.Losses.mse(model(x), y)

# 5. Define an optimizer
#    Simple Gradient Descent with a learning rate of 0.1
opt = Flux.Optimise.Descent(0.1)

# 6. Train the model for a number of epochs
epochs = 200
println("\nStarting training...")
for epoch in 1:epochs
    Flux.train!(loss, Flux.params(model), data, opt)
    # Optional: Print loss every few epochs
    if epoch % 50 == 0
        current_loss = loss(xtrain, ytrain)
        println("Epoch: $epoch, Loss: $current_loss")
    end
end
println("Training finished.")

# 7. Check the learned parameters
#    They should be close to W=[2.0] and b=[-1.0]
println("\nLearned parameters: ", Flux.params(model))

# 8. Make a prediction
x_test = hcat([5.0f0]) # Test with a new input (x=5)
prediction = model(x_test)
println("\nPrediction for x=5: ", prediction)
println("(Expected value is 2*5 - 1 = 9)") 