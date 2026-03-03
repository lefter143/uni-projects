prev_wumpa_fruit = 0
prev_crates = 0
prev_lives = 0
prev_crash_x = 0
prev_enemies_alive = 0

-- Enemy killing confirmation timer (frames)(instead of rendering and derendering)
frames_to_kill_counter = 0
frames_to_kill = 90  -- 1.5 seconds at ~60 FPS


moving_factor = 0.03
time_penalty = 0.001

-- Calculate rewards
function total_reward()
    reward_wumpa_fruit = wumpa_fruit_reward()
    reward_crates = crates_reward()
    reward_from_enemy = enemy_reward()
    penalty_lives = lives_penalty()
    reward_for_moving = moving_reward() * moving_factor
    reward_time_penalty = -time_penalty

    -- Clamp moving reward to reasonable range to avoid large spikes
    if reward_for_moving < -1 then 
        reward_for_moving = -moving_factor
    elseif reward_for_moving > 1 then
        reward_for_moving = moving_factor
    end

    -- Handle edge cases where multiple rewards might explode
    if reward_crates < 0 and reward_wumpa_fruit < 0 and penalty_lives < -1 then
        reward_crates = 0
        reward_wumpa_fruit = 0
        penalty_lives = 0
    end
    if reward_crates + reward_wumpa_fruit + penalty_lives > 10 then
        reward_crates = 0
        reward_wumpa_fruit = 0
        penalty_lives = 0
    end

    -- Update previous states
    prev_wumpa_fruit = data.wumpa_fruit 
    prev_crates = data.crates
    prev_lives = data.lives
    prev_crash_x = data.crash_x

    return reward_wumpa_fruit + reward_crates + reward_from_enemy + penalty_lives + reward_for_moving + reward_time_penalty
end

function wumpa_fruit_reward()
    diff = data.wumpa_fruit - prev_wumpa_fruit
    if diff < 0 then
        -- Handle mod 100(0 after 99)
        diff = 1
    end
    prev_wumpa_fruit = data.wumpa_fruit
    return diff
end

function crates_reward()
    diff = data.crates - prev_crates
    if diff < 0 then
	prev_crates = data.crates
	diff=0
    end
    return diff
end


function lives_penalty()
    return 2*(data.lives - prev_lives)
end

function moving_reward()
    return data.crash_x - prev_crash_x
end

function enemy_reward()
    current = data.enemies_tokill  -- Assuming this is your variable counting alive enemies

    -- If enemy count decreased, reset confirmation timer
    if current < prev_enemies_alive then
        frames_to_kill_counter = frames_to_kill
        prev_enemies_alive = current
    elseif frames_to_kill_counter > 0 then
        frames_to_kill_counter = frames_to_kill_counter - 1
    end

    -- Only confirm enemy kill reward when timer expires and enemy count stayed low
    if frames_to_kill_counter == 0 and current < prev_enemies_alive then
        prev_enemies_alive = current
        return 1
    end

    return 0
end
